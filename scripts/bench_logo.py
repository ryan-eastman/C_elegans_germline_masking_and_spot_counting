"""Leave-one-gonad-out estimates of the spot detectors (docs/BENCHMARK_PROTOCOL.md, section 6.2;
docs/BENCHMARK.md, "Leave-one-gonad-out"). Train and val gonads only; a test gonad is never opened.

    python scripts/bench_logo.py folds
    python scripts/bench_logo.py select --backend log|spotmax --runs DIR ... --out FILE [--search d3|prereg-grid]
    python scripts/bench_logo.py grid-configs --out FILE
    python scripts/bench_logo.py heldout --detector NAME --backend B (--params JSON | --selection FILE)
        --out-root DIR [--tags TAG ...]
    python scripts/bench_logo.py net3d-fold --held-out TAG --recipe-from MODELDIR --models-root DIR
        --out-root DIR [--device cuda] [--prep-dir DIR] [--min-free-gb 8]
    python scripts/bench_logo.py net3d-folds --recipe-from MODELDIR --models-root DIR --out-root DIR
        --log-dir DIR [--gpu-lock LOCK.py] [--tags TAG ...] (the same options as net3d-fold)
    python scripts/bench_logo.py report --heldout NAME=DIR ... --out DIR [--tags TAG ...]

folds        the six folds of section 6.2 (germquant.benchmark.logo.logo_folds) as JSON.
select       replays the D3 search of a tuned classical detector on every fold's own gonads
             (germquant.benchmark.logo.replay_log_search / replay_spotmax_search) from the bench_tune.py
             and bench_detector.py run folders given with --runs, and writes the selected parameter set of
             every fold, or the parameter sets and gonads still to be run (exit code 3).
             --search prereg-grid (backend spotmax; --search d3, the D3 replay, is the default) applies the
             rule of the pre-registered SpotMAX grid instead (protocol section 8 item 2, sections 6.1 and
             6.4; germquant.benchmark.logo.replay_grid). On every fold and on the designated split (the
             split's own train and val gonads), all 60 sets are ranked on the selection gonads by the
             primary metric (the truth-mask CCC rounded to 3 decimals, then the spot F1 inside R, then the
             grid's order) and the best is selected. The file has the layout of a D3 selection, which
             heldout --selection reads, plus "search", "configs" (the grid), "designated" (the designated
             split's selection) and "to_run" (the missing pairs merged by set). Exit code 3 as above; a
             score without the truth-mask counts is missing too. Run folders that disagree, or a fold
             without selection gonads, are refused (exit code 2).
grid-configs writes the 60 sets of the pre-registered SpotMAX grid
             (germquant.benchmark.logo.prereg_spotmax_grid) as a configuration file for
             scripts/bench_tune.py --configs.
heldout      scores the held-out gonad of every fold (or of --tags) with scripts/bench_detector.py
             --detector-only: one run, one trial record (with a "note" naming the fold) and one output folder
             <out-root>/<tag>/ per gonad. --params gives one parameter set for every fold (the deployed
             SpotMAX recipe and the log defaults are {}); --selection the per-fold sets of `select`.
net3d-fold   trains net3d on the fold's five gonads by the frozen recipe of a selected model (germquant
             train-spots --recipe-from: its settings, seed, schedule, epoch count and threshold, no
             validation) into <models-root>/fold_<held-out>/, then scores the held-out gonad as `heldout`
             does. Run it under the workstation's GPU lock.
net3d-folds  runs net3d-fold for every fold (or --tags) in turn, each in its own process and, with
             --gpu-lock, each under that lock script; logs in --log-dir.
report       pools the held-out gonads of every detector (cluster bootstrap by gonad, 1000 resamples, seed
             20251105) and compares every pair of detectors on the same resamples; writes logo_report.json,
             logo_per_gonad.csv, logo_pooled.csv and logo_paired.csv. With --tags it pools those gonads from
             the folders given (each DIR may then be one bench_detector.py output folder), which is how the
             one-time test scoring of docs/FINAL_SCORING.md compares the competitors; it computes nothing
             itself beyond pooling the per-gonad and per-nucleus tables it is given.

Common options: --manifest, --cache-root, --trials-log, --split-file (net3d-fold), --n-boot, --seed.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd  # noqa: E402

from germquant.benchmark import detector as bd  # noqa: E402
from germquant.benchmark import logo  # noqa: E402
from germquant.benchmark.manifest import load_manifest  # noqa: E402
from germquant.benchmark.report import jsonable  # noqa: E402

DEFAULT_SPLIT_FILE = REPO / "tests" / "golden" / "rad51_split.json"


def _bench_detector():
    spec = importlib.util.spec_from_file_location("bench_detector", REPO / "scripts" / "bench_detector.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _folds(args) -> list[dict]:
    return logo.logo_folds(load_manifest(args.manifest), cache_root=args.cache_root)


def _fold_of(folds: list[dict], tag: str) -> dict:
    for f in folds:
        if f["held_out"] == tag:
            return f
    raise SystemExit(f"REFUSED: {tag} is not a held-out gonad of any fold (train or val gonads with usable spot "
                     f"truth: {[f['held_out'] for f in folds]})")


def cmd_folds(args) -> int:
    print(json.dumps(_folds(args), indent=1))
    return 0


def cmd_select(args) -> int:
    if args.search == "prereg-grid":
        return select_prereg_grid(args)
    idx = logo.index_runs(args.runs, backend=args.backend)
    replay = logo.replay_log_search if args.backend == "log" else logo.replay_spotmax_search
    folds = _folds(args)
    out = {"backend": args.backend, "runs": [str(r) for r in args.runs], "objective": logo.bt.OBJECTIVE_TEXT,
           "folds": []}
    missing = []
    for f in folds:
        r = replay(idx, f["train"], f["val"])
        out["folds"].append({**f, **r})
        if not r.get("complete"):
            missing += [{**m, "fold": f["held_out"], "stage": r.get("stage")} for m in r["missing"]]
        print(f"{f['held_out']}: train {f['train']}, val {f['val']} -> "
              + (f"selected {r['selected']} (finalists {r['finalists']})" if r.get("complete")
                 else f"incomplete at stage {r.get('stage')}: {len(r['missing'])} sets to run"))
    out["complete"] = not missing
    out["missing"] = missing
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(jsonable(out), indent=1) + "\n", encoding="utf-8")
    print(f"-> {args.out}")
    return 0 if not missing else 3


def select_prereg_grid(args) -> int:
    """``select --search prereg-grid``: the rule of the pre-registered SpotMAX grid (protocol sections 6.1,
    6.4 and 8 item 2) on every fold and on the designated split, written in the layout of the D3 replay.
    Run folders that disagree, a fold without selection gonads, and a fold where no set has a finite
    truth-mask CCC are refused (exit code 2); nothing is written then."""
    if args.backend != "spotmax":
        print("REFUSED: --search prereg-grid is the SpotMAX grid of protocol section 8 item 2 (--backend spotmax)")
        return 2
    configs = logo.prereg_spotmax_grid()
    out = {"backend": args.backend, "search": args.search, "runs": [str(r) for r in args.runs],
           "objective": logo.PRIMARY_OBJECTIVE_TEXT, "configs": configs, "folds": []}
    missing = []

    def replay(label: str, unit: dict) -> dict:
        r = logo.replay_grid(idx, configs, unit["train"], unit["val"])
        if r["complete"]:
            ccc, f1 = r["val_objective"][r["selected"]]
            print(f"{label}: val {unit['val']} -> selected {r['selected']} (val truth-mask CCC {ccc:.3f}, "
                  f"F1 inside R {f1:.3f}; {r['n_sets']} sets)")
        else:
            print(f"{label}: val {unit['val']} -> incomplete: {len(r['missing'])} sets lack a score with the "
                  "truth-mask counts on a val gonad")
            missing.extend({**m, "fold": label, "stage": r["stage"]} for m in r["missing"])
        return r

    try:
        folds = _folds(args)
        idx = logo.index_runs(args.runs, backend="spotmax", truth_masks=True)
        for f in folds:
            out["folds"].append({**f, **replay(f["held_out"], f)})
        designated = logo.designated_split(folds)
        out["designated"] = {**designated, **replay("designated", designated)}
    except ValueError as err:
        print(f"REFUSED: {err}")
        return 2
    out["complete"] = not missing
    out["missing"] = missing
    to_run: dict[str, set] = {}
    for m in missing:
        to_run.setdefault(m["name"], set()).update(m["tags"])
    out["to_run"] = [{"name": c["name"], "params": c["params"], "tags": sorted(to_run[c["name"]])}
                     for c in configs if c["name"] in to_run]
    if out["to_run"]:
        per_tag: dict[str, int] = {}
        for m in out["to_run"]:
            for t in m["tags"]:
                per_tag[t] = per_tag.get(t, 0) + 1
        print("still to run: " + ", ".join(f"{t} ({n} sets)" for t, n in sorted(per_tag.items())))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(jsonable(out), indent=1) + "\n", encoding="utf-8")
    print(f"-> {args.out}")
    return 0 if not missing else 3


def cmd_grid_configs(args) -> int:
    configs = logo.prereg_spotmax_grid()
    doc = {"backend": "spotmax", "search": logo.PREREG_SPOTMAX_SEARCH,
           "source": "germquant.benchmark.logo.prereg_spotmax_grid, docs/BENCHMARK_PROTOCOL.md section 8 item 2",
           "configs": configs}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(f"{len(configs)} parameter sets of the pre-registered SpotMAX grid -> {out}")
    return 0


def _common_bench_args(args) -> list[str]:
    a = ["--manifest", str(args.manifest), "--cache-root", str(args.cache_root), "--n-boot", str(args.n_boot),
         "--seed", str(args.seed)]
    if args.trials_log:
        a += ["--trials-log", str(args.trials_log)]
    return a


def score_heldout(args, fold: dict, detector: str, backend: str, params: dict, extra_note: dict | None = None) -> int:
    out = Path(args.out_root) / fold["held_out"]
    note = {"logo": {"detector": detector, "held_out": fold["held_out"], "fold_train": fold["train"],
                     "fold_val": fold["val"], "substitute": fold["substitute"], **(extra_note or {})}}
    argv = ["--backend", backend, "--params", json.dumps(params), "--split", fold["role"], "--tags", fold["held_out"],
            "--detector-only", "--out", str(out), "--note", json.dumps(note), *_common_bench_args(args)]
    print(f"[{detector}] held out {fold['held_out']}: params {bd.canonical_json(params)}", flush=True)
    return _bench_detector().main(argv)


def cmd_heldout(args) -> int:
    folds = _folds(args)
    if args.tags:
        folds = [_fold_of(folds, t) for t in args.tags]
    selection = None
    if args.selection:
        selection = json.loads(Path(args.selection).read_text(encoding="utf-8"))
        if not selection.get("complete"):
            print("REFUSED: the selection file is incomplete (run its missing parameter sets first)")
            return 2
        by_tag = {f["held_out"]: f for f in selection["folds"]}
    elif args.params is None:
        print("REFUSED: give --params or --selection")
        return 2
    rc = 0
    for f in folds:
        if selection is not None:
            s = by_tag[f["held_out"]]
            if s["train"] != f["train"] or s["val"] != f["val"]:
                print(f"REFUSED: the selection file's fold for {f['held_out']} has other gonads")
                return 2
            params, extra = s["selected_params"], {"selected": s["selected"], "finalists": s["finalists"],
                                                   "selection_file": str(args.selection)}
        else:
            params, extra = json.loads(args.params), {"params_for_every_fold": True}
        rc = max(rc, score_heldout(args, f, args.detector, args.backend, params, extra))
    return rc


def cmd_net3d_fold(args) -> int:
    from germquant.spots.net3d.train import recipe_from, train_spots

    folds = _folds(args)
    fold = _fold_of(folds, args.held_out)
    kw, recipe = recipe_from(args.recipe_from)
    # the five gonads, in the recipe's order: its train gonads, then its val gonads, the held-out one left out
    order = [t for t in recipe["from_train_gonads"] + recipe["from_val_gonads"] if t != fold["held_out"]]
    fold_tags = set(fold["train"]) | set(fold["val"])
    if set(order) != fold_tags:
        print(f"REFUSED: the recipe's gonads {order} are not the fold's {sorted(fold_tags)}")
        return 2
    model_dir = Path(args.models_root) / f"fold_{fold['held_out']}"
    if (model_dir / "model.pt").is_file() and not args.retrain:
        print(f"{model_dir / 'model.pt'} exists; scoring it (pass --retrain to train again)")
    else:
        res = train_spots(args.cache_root, args.split_file, model_dir, manifest=args.manifest, device=args.device,
                          prep_dir=args.prep_dir, min_free_gb=args.min_free_gb, train_tags=order, val_tags=[],
                          recipe=recipe, name=f"net3d_logo_fold_{fold['held_out']}", **kw)
        print(f"trained {res['model_path']} (sha256 {res['model_sha256'][:12]}...) on {order}: {res['stop_reason']}")
    from germquant.spots.net3d.infer import file_sha256

    model = (model_dir / "model.pt").as_posix()
    params = {"net3d": {"model": model, "device": args.device}}
    return score_heldout(args, fold, "net3d", "net3d", params,
                         {"model": model, "model_sha256": file_sha256(model), "trained_on": order,
                          "recipe_from": recipe["from_model_dir"], "recipe_model_sha256": recipe["from_model_sha256"]})


def cmd_net3d_folds(args) -> int:
    """Every fold's ``net3d-fold`` in turn, each in its own process; with --gpu-lock each one runs under
    that lock script (``LOCK run --name NAME -- COMMAND``), so other jobs can take the GPU between folds."""
    import subprocess

    folds = _folds(args)
    tags = list(args.tags) if args.tags else [f["held_out"] for f in folds]
    for t in tags:
        _fold_of(folds, t)
    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    common = ["--manifest", str(args.manifest), "--cache-root", str(args.cache_root), "--n-boot", str(args.n_boot),
              "--seed", str(args.seed)] + (["--trials-log", str(args.trials_log)] if args.trials_log else [])
    rc = 0
    for t in tags:
        cmd = [sys.executable, str(Path(__file__).resolve()), *common, "net3d-fold", "--held-out", t,
               "--recipe-from", str(args.recipe_from), "--models-root", str(args.models_root), "--out-root",
               str(args.out_root), "--split-file", str(args.split_file), "--device", args.device,
               "--min-free-gb", str(args.min_free_gb)] + (["--prep-dir", str(args.prep_dir)] if args.prep_dir else []) \
            + (["--retrain"] if args.retrain else [])
        if args.gpu_lock:
            cmd = [sys.executable, str(args.gpu_lock), "run", "--name", f"net3d LOGO fold held-out {t}", "--", *cmd]
        with open(log_dir / f"net3d_fold_{t}.out.log", "w", encoding="utf-8") as fo, \
                open(log_dir / f"net3d_fold_{t}.err.log", "w", encoding="utf-8") as fe:
            code = subprocess.run(cmd, stdout=fo, stderr=fe).returncode
        print(f"fold {t}: exit {code}", flush=True)
        rc = max(rc, code)
    return rc


def cmd_report(args) -> int:
    if args.tags:
        folds, tags = [], list(args.tags)
    else:
        folds = _folds(args)
        tags = [f["held_out"] for f in folds]
    held = {}
    for spec in args.heldout:
        name, _, root = spec.partition("=")
        root = Path(root)
        if (root / "per_gonad.csv").is_file():
            dirs = [root]                                   # one bench_detector.py folder holding every gonad
        else:
            dirs = [root / t for t in tags]                 # one folder per held-out gonad
        lack = [str(d) for d in dirs if not (d / "per_gonad.csv").is_file()]
        if lack:
            print(f"REFUSED: {name} has no held-out scores in {lack}")
            return 2
        held[name] = logo.load_heldout(name, dirs)
        have = set(held[name].per_gonad["tag"].astype(str))
        if not set(tags) <= have:
            print(f"REFUSED: {name} has no scores for {sorted(set(tags) - have)}")
            return 2
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    pooled = {n: logo.pool_heldout(h, tags, n_boot=args.n_boot, seed=args.seed) for n, h in held.items()}
    pairs = {f"{a} - {b}": logo.paired(held[a], held[b], tags, n_boot=args.n_boot, seed=args.seed)
             for a, b in itertools.combinations(list(held), 2)}
    per_gonad = pd.concat([logo.per_gonad_table(h).assign(detector=n) for n, h in held.items()], ignore_index=True)
    per_gonad.to_csv(out / "logo_per_gonad.csv", index=False)
    rows = []
    for n, blocks in pooled.items():
        for block, rec in blocks.items():
            for stat, s in rec["stats"].items():
                rows.append({"detector": n, "block": block, "stat": stat, "point": s["point"], "lo": s["lo"],
                             "hi": s["hi"], "n_nuclei": rec.get("n_nuclei"), "n_gonads": rec["n_gonads"]})
    pd.DataFrame(rows).to_csv(out / "logo_pooled.csv", index=False)
    prow = []
    for pair, blocks in pairs.items():
        for block, stats_ in blocks.items():
            for stat, s in stats_.items():
                prow.append({"pair": pair, "block": block, "stat": stat, "a": s["a"], "b": s["b"],
                             "diff": s["diff"]["point"], "lo": s["diff"]["lo"], "hi": s["diff"]["hi"],
                             "excludes_zero": s["excludes_zero"]})
    pd.DataFrame(prow).to_csv(out / "logo_paired.csv", index=False)
    from germquant import provenance

    doc = {"folds": folds, "held_out_gonads": tags, "detectors": {n: h.sources for n, h in held.items()},
           "pooled": pooled, "paired": pairs, "bootstrap": {"n_boot": args.n_boot, "seed": args.seed, "unit": "gonad"},
           "git_sha": provenance.git_sha()}
    (out / "logo_report.json").write_text(json.dumps(jsonable(doc), indent=1) + "\n", encoding="utf-8")
    for n, blocks in pooled.items():
        tm, g, f1 = blocks.get("truth_masks"), blocks.get("germline"), blocks.get("spots_all")
        f1r = blocks.get("spots_truth_masks")

        def fmt(rec, stat):
            if not rec:
                return "n/a"
            s = rec["stats"][stat]
            return f"{s['point']:.3f} [{s['lo']:.3f}, {s['hi']:.3f}]"

        print(f"{n}: truth masks CCC {fmt(tm, 'ccc')} | germline CCC {fmt(g, 'ccc')} | F1 {fmt(f1, 'f1')} | "
              f"F1 inside R {fmt(f1r, 'f1')}")
    print(f"-> {out}")
    return 0


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=str(bd.DEFAULT_MANIFEST))
    ap.add_argument("--cache-root", default=str(bd.DEFAULT_CACHE_ROOT))
    ap.add_argument("--trials-log", default=None, help="default: trials.jsonl next to the manifest")
    ap.add_argument("--n-boot", type=int, default=bd.N_BOOT)
    ap.add_argument("--seed", type=int, default=bd.BOOT_SEED)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("folds")
    s = sub.add_parser("select")
    s.add_argument("--backend", required=True, choices=["log", "spotmax"])
    s.add_argument("--runs", nargs="+", required=True, help="bench_tune.py / bench_detector.py output folders")
    s.add_argument("--out", required=True)
    s.add_argument("--search", default="d3", choices=["d3", "prereg-grid"],
                   help="d3: replay the D3 search (default); prereg-grid: the pre-registered SpotMAX grid "
                        "(protocol section 8 item 2)")
    g = sub.add_parser("grid-configs")
    g.add_argument("--out", required=True, help="JSON file for scripts/bench_tune.py --configs")
    h = sub.add_parser("heldout")
    h.add_argument("--detector", required=True, help="a name for the detector (the report's row label)")
    h.add_argument("--backend", required=True)
    h.add_argument("--params", default=None, help="one JSON parameter set for every fold")
    h.add_argument("--selection", default=None, help="the per-fold selection file of `select`")
    h.add_argument("--out-root", required=True)
    h.add_argument("--tags", nargs="*", default=None, help="only these held-out gonads")
    n = sub.add_parser("net3d-fold")
    n.add_argument("--held-out", required=True)
    n.add_argument("--recipe-from", required=True)
    n.add_argument("--models-root", required=True)
    n.add_argument("--out-root", required=True)
    n.add_argument("--split-file", default=str(DEFAULT_SPLIT_FILE))
    n.add_argument("--device", default="cuda", choices=["cuda", "cpu", "auto"])
    n.add_argument("--prep-dir", default=None)
    n.add_argument("--min-free-gb", type=float, default=0.0)
    n.add_argument("--retrain", action="store_true")
    na = sub.add_parser("net3d-folds")
    na.add_argument("--recipe-from", required=True)
    na.add_argument("--models-root", required=True)
    na.add_argument("--out-root", required=True)
    na.add_argument("--log-dir", required=True)
    na.add_argument("--tags", nargs="*", default=None, help="only these held-out gonads (default: every fold)")
    na.add_argument("--gpu-lock", default=None, help="a lock script run as LOCK run --name NAME -- COMMAND")
    na.add_argument("--split-file", default=str(DEFAULT_SPLIT_FILE))
    na.add_argument("--device", default="cuda", choices=["cuda", "cpu", "auto"])
    na.add_argument("--prep-dir", default=None)
    na.add_argument("--min-free-gb", type=float, default=0.0)
    na.add_argument("--retrain", action="store_true")
    r = sub.add_parser("report")
    r.add_argument("--heldout", nargs="+", required=True,
                   help="NAME=DIR, DIR holding <tag>/ held-out folders or one bench_detector.py output folder")
    r.add_argument("--tags", nargs="*", default=None,
                   help="pool these gonads instead of the folds' held-out gonads (every DIR must hold them all)")
    r.add_argument("--out", required=True)
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    os.environ.setdefault("TQDM_DISABLE", "1")
    return {"folds": cmd_folds, "select": cmd_select, "grid-configs": cmd_grid_configs, "heldout": cmd_heldout,
            "net3d-fold": cmd_net3d_fold, "net3d-folds": cmd_net3d_folds, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
