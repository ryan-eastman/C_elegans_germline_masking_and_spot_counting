"""Spot preset tuning: search a grid of spot keys on finished runs and rank the parameter sets by a fixed
objective against truth points (design 6.8 and 10 step 9; plan D24; docs/TUNING.md).

Spots are tuned, not trained. A search takes finished runs of annotated images (their label images, their
geometry and the config they were made with) and runs only the spot detection again, once per parameter set
of the grid (``germquant.tune.grid``), with the same call the spots stage of process_image makes: every key
the grid does not name keeps the value of the config the run was made with. Each result is scored against the
image's truth points (``germquant.tune.points``) inside the annotated objects, the objects that hold at least
one truth point:

* per object, the truth count (its points) against the detected count, pooled over the objects of every
  image: Lin's CCC (``object_ccc``). Two equal series of at least two objects count as a CCC of 1.0: the
  formula is 0 / 0 when every annotated object holds the same count, and the set that finds exactly those
  counts must not rank below every other set (``_ccc``);
* per point, a one-to-one matching of the detected spots inside those objects with the truth points inside
  an anisotropic cylinder of 0.3 um in xy and 0.4 um in z (the RAD-51 harness's tolerances; points placed on
  a projection are matched in 2D): precision, recall and F1, pooled from the summed counts (``point_f1``).

Objects without points are not scored, so a detection in an object nobody annotated is neither a false
positive nor a miss, and the search tunes the detector, not the segmentation (the labels stay fixed). The
objective is the pooled object CCC rounded to 3 decimals, then the pooled point F1 (higher is better; a
missing value ranks last; equal keys keep the grid's order): the order the RAD-51 harness ranks by (D24).

Two tuners share that code:

* ``tune_spotmax`` runs SpotMAX once per group of parameter sets that differ only in
  ``spots.effect_size_min``, with the minimum at 0, and applies each set's minimum afterwards with
  ``germquant.tune.effect_filter.filter_effect_size`` (design 6.8). sharpen is part of the group key like
  every other key, so a threshold is never carried from one sharpen value to another;
* ``tune_grid`` runs any backend (log included) once per parameter set.

A tuning run writes its folder (``out_dir``; tuning/<run_id>/ in a project): grid.yaml (the grid searched,
in the --grid format, so the search can be repeated), trials.jsonl (one line per parameter set, with the
per-image and pooled scores and the point and object counts), ranked.csv (rank, objective, objective_value,
tie_break, tie_break_value, then one column per grid key; values to 3 decimals), best.yaml ({objective,
backend, keys}; not written when no set could be scored) and run.json (what the run read: the finished runs,
their config hashes and done-marker digests, the point files and their digests, taken on LF-normalised bytes,
the tolerances, the engine version). Every file is ASCII with LF line endings, written atomically.

Imports (plan D28): the module imports nothing heavy. Every tuner imports germquant.pipeline (torch) before
it reads an image or loads a backend, so torch claims cuBLAS before SpotMAX or CuPy can (``import spotmax.pipe``
alone loads CuPy without torch). Of the benchmark and validation code only germquant.benchmark.stats,
matching, nuclei3d and spots and germquant.validate.compare may be imported (the seal of design 6.7;
tests/test_tune_seal.py).
"""
from __future__ import annotations

import copy
import json
import logging
import math
import os
import re
import secrets
import time
from pathlib import Path

import numpy as np
import pandas as pd

from germquant.fsutil import long_path

from .effect_filter import filter_effect_size
from .grid import (
    EFFECT_SIZE_KEY,
    TuneRefused,
    check_grid,
    default_grid,
    group_key,
    load_grid,
    n_groups,
    n_runs,
    parameter_sets,
)
from .points import check_points, parse_points_name, read_points, to_run_frame, write_text_atomic

log = logging.getLogger(__name__)

CCC_DECIMALS = 3
METRICS = ("object_ccc", "point_f1", "point_precision", "point_recall")
OBJECTIVE = ("object_ccc", "point_f1")
_METRIC_TEXT = {"object_ccc": "per-object Lin's CCC over the annotated objects",
                "point_f1": "point F1", "point_precision": "point precision", "point_recall": "point recall"}
TOL_XY_UM = 0.3              # the RAD-51 harness's matching cylinder (germquant.benchmark.spots)
TOL_Z_UM = 0.4
MAX_RUNS = 50
FILES = {"grid": "grid.yaml", "trials": "trials.jsonl", "ranked": "ranked.csv", "best": "best.yaml",
         "run": "run.json"}
RUN_ID_RE = re.compile(r"\d{8}T\d{6}Z_[0-9a-f]{8}")
_SPLIT_RE = re.compile(r"[A-Za-z0-9_-]{1,32}")


def objective_text(objective=OBJECTIVE) -> str:
    """The objective as it is recorded in best.yaml and run.json."""
    primary, tie = _check_objective(objective)
    return (f"{_METRIC_TEXT[primary]}, rounded to {CCC_DECIMALS} decimals, then {_METRIC_TEXT[tie]} "
            f"(higher is better; a missing value ranks last)")


def _check_objective(objective) -> tuple[str, str]:
    try:
        primary, tie = objective
    except (TypeError, ValueError):
        raise TuneRefused(f"objective must be a pair of metric names from {METRICS}, got {objective!r}") from None
    if primary not in METRICS or tie not in METRICS or primary == tie:
        raise TuneRefused(f"objective must be two different metrics from {METRICS}, got {objective!r}")
    return str(primary), str(tie)


OBJECTIVE_TEXT = objective_text(OBJECTIVE)


def new_run_id() -> str:
    """A tuning run id, <UTC time>_<8 hex>, in the form of the run-log ids (^\\d{8}T\\d{6}Z_[0-9a-f]{8}$)."""
    return time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()) + "_" + secrets.token_hex(4)


# ------------------------------------------------------------------------------ public tuners
def tune_spotmax(items, grid, objective=OBJECTIVE, split="project", *, out_dir, max_runs=MAX_RUNS) -> dict:
    """Tune the SpotMAX keys of ``grid`` on ``items`` (finished runs with truth points), running SpotMAX once
    per group of parameter sets that differ only in ``spots.effect_size_min`` and filtering each set's
    minimum afterwards (design 6.8).

    ``items`` is a list of {"image_id", "image" (the image file), "result_dir" (the image's finished run:
    labels, done marker, run_manifest.json, image summary), "points" (a truth-point file or the DataFrame
    ``read_points`` returns)}; ``grid`` a grid ({dotted spots key: [values]}); ``objective`` a pair of metric
    names (default: object CCC, then point F1); ``split`` the name recorded with the run (a benchmark split,
    or "project"; "test" is refused: tuning never reads a test split); ``out_dir`` the run folder to write
    (absent or empty); ``max_runs`` the most SpotMAX runs allowed. Returns {"run_dir", "best" (best.yaml's
    content, or None when no set could be scored), "ranked" (ranked.csv's rows), "n_runs", "warnings"}.
    Raises TuneRefused before any detection for a grid, an item or a budget it cannot run."""
    return _tune(items, grid, objective, split, backend="spotmax", out_dir=out_dir, max_runs=max_runs,
                 shortcut=True)


def tune_grid(items, grid, objective=OBJECTIVE, split="project", *, backend, out_dir, max_runs=MAX_RUNS) -> dict:
    """Tune the keys of ``grid`` for the spot backend ``backend`` (log, spotmax or one from an entry point) on
    ``items``, with one detection per parameter set and image (D24). The arguments and the result are those of
    ``tune_spotmax``; ``max_runs`` bounds the detections."""
    return _tune(items, grid, objective, split, backend=str(backend), out_dir=out_dir, max_runs=max_runs,
                 shortcut=False)


def tune_spots(project, points="annotations", grid=None, *, out=None, max_runs=MAX_RUNS) -> dict:
    """Tune the spot keys of a project's recipe on the truth points of its annotated images (design 10 step 9;
    ``germquant tune-spots``).

    ``project`` is a germquant.project.Project (or its folder); ``points`` the folder of truth-point files,
    relative to the project (default annotations/), holding one ``<image_id>__points_<annotator>.csv`` per
    image (two annotators of one image are refused: keep each in its own folder); ``grid`` a grid, a grid file
    (relative to the project) or None for the default grid of the recipe's spot backend; ``out`` the run
    folder (default tuning/<run_id>/, relative to the project). Each annotated image needs a finished run in
    results/<image_id>/. The recipe's backend is tuned: spotmax through ``tune_spotmax``, any other through
    ``tune_grid``. Returns {"run_dir", "best", "ranked", "warnings"}; apply_tuning writes the best set into
    recipe.yaml."""
    import germquant.pipeline  # noqa: F401  (D28: torch before any reader, backend or CuPy)

    open_project, project_config = _project_api()
    if isinstance(project, (str, os.PathLike)):
        project = open_project(project)
    try:
        cfg = project_config(project)
    except (ValueError, FileNotFoundError) as e:        # germquant.project's ProjectError is a ValueError
        raise TuneRefused(f"the project's recipe cannot be loaded: {e}") from e
    root = Path(project.root)
    backend = str(cfg.get("spots.backend", "spotmax"))
    if grid is None:
        grid = default_grid(cfg)
    elif isinstance(grid, (str, os.PathLike)):
        gp = Path(os.fspath(grid))
        grid = load_grid(gp if gp.is_absolute() else root / gp)
    items, warnings = project_items(project, points)
    run_id = new_run_id()
    out_dir = root / "tuning" / run_id if out is None else Path(os.fspath(out))
    if not out_dir.is_absolute():
        out_dir = root / out_dir
    current = {"config_hash": cfg.hash, "spots": (cfg.as_dict().get("spots") or {})}
    res = _tune(items, grid, OBJECTIVE, "project", backend=backend, out_dir=out_dir, max_runs=max_runs,
                shortcut=backend == "spotmax", warnings=warnings, current=current, run_id=run_id, root=root)
    return {"run_dir": res["run_dir"], "best": res["best"], "ranked": res["ranked"], "warnings": res["warnings"]}


def _project_api():
    """open_project and project_config of germquant.project (a seam the tests replace)."""
    from germquant import project as P

    return P.open_project, P.project_config


def project_items(project, points="annotations") -> tuple[list[dict], list[str]]:
    """The tuning items of a project: one per image of its manifest that has a truth-point file in the
    ``points`` folder (relative to the project), with its finished run in results/<image_id>/. Returns
    (items, warnings). Files whose names are not truth-point names, files of images the project does not
    hold, and excluded or time-series images are skipped with a warning; two files for one image (two
    annotators) and a folder without any usable file are refused (TuneRefused)."""
    root = Path(project.root)
    pdir = Path(os.fspath(points))
    pdir = pdir if pdir.is_absolute() else root / pdir
    if not os.path.isdir(long_path(pdir)):
        raise TuneRefused(f"no annotations folder {pdir}: save truth points as <image_id>__points_<annotator>.csv "
                          f"there first (docs/TUNING.md)")
    warnings: list[str] = []
    by_image: dict[str, list[Path]] = {}
    names = sorted(n for n in os.listdir(long_path(pdir)) if n.lower().endswith(".csv"))
    for f in (pdir / n for n in names if os.path.isfile(long_path(pdir / n))):
        parsed = parse_points_name(f.name)
        if parsed is None:
            warnings.append(f"{f.name}: not a truth-point file name (<image_id>__points_<annotator>.csv); skipped")
            continue
        by_image.setdefault(parsed[0], []).append(f)
    images = {str(im.get("id")): im for im in (project.images or [])}
    items = []
    for image_id in sorted(by_image):
        files = by_image[image_id]
        if len(files) > 1:
            raise TuneRefused(f"{image_id} has truth points of several annotators ({', '.join(f.name for f in files)}); "
                              f"tuning scores against one: keep each annotator's files in a folder of their own "
                              f"and pass it with --points")
        im = images.get(image_id)
        if im is None:
            warnings.append(f"{files[0].name}: {image_id} is not an image of this project; skipped")
            continue
        if im.get("excluded"):
            warnings.append(f"{files[0].name}: {image_id} is excluded from the project; skipped")
            continue
        if im.get("timepoints"):
            warnings.append(f"{files[0].name}: {image_id} is a time series, which spot tuning does not take in this "
                            f"release; skipped")
            continue
        path = Path(str(im.get("path")))
        items.append({"image_id": image_id, "image": path if path.is_absolute() else root / path,
                      "result_dir": root / "results" / image_id, "points": files[0]})
    if not items:
        raise TuneRefused(f"no truth-point file in {pdir} belongs to an image of this project"
                          + (f" ({'; '.join(warnings)})" if warnings else ""))
    return items, warnings


# ------------------------------------------------------------------------------ the search
def _tune(items, grid, objective, split, *, backend, out_dir, max_runs, shortcut, warnings=None, current=None,
          run_id=None, root=None, on_tables=None) -> dict:
    """The search behind both tuners. ``shortcut`` runs SpotMAX once per effect-size group; ``current`` (the
    project's config hash and spots: section) adds a warning for each run made with another config;
    ``on_tables(set_index, item_index, per_spot, per_nucleus)`` sees every detection's tables (tests)."""
    import germquant.pipeline  # noqa: F401  (D28: torch before any reader, backend or CuPy)

    t_start = time.time()
    primary, tie = _check_objective(objective)
    split = str(split)
    if not _SPLIT_RE.fullmatch(split):
        raise TuneRefused(f"split {split!r}: use 1 to 32 letters, digits, '-' or '_'")
    if split.lower() == "test":
        raise TuneRefused("tuning never runs on a test split: tune on the tune split, select on val "
                          "(design 6.6)")
    grid = check_grid(grid, backend if backend in ("spotmax", "log") else None)
    _check_backend(backend)
    if not isinstance(items, (list, tuple)) or not items:
        raise TuneRefused("tuning needs at least one item (a finished run with truth points)")
    max_runs = int(max_runs)
    total = n_runs(grid, len(items), backend if shortcut else "")
    if total > max_runs:
        per = n_groups(grid, backend if shortcut else "")
        raise TuneRefused(f"the grid means {total} {backend} runs ({per} per image on {len(items)} image(s)), above "
                          f"--max-runs {max_runs}; narrow the grid or raise --max-runs")
    out_dir = Path(os.fspath(out_dir))
    if os.path.exists(long_path(out_dir)) and (not os.path.isdir(long_path(out_dir))
                                               or os.listdir(long_path(out_dir))):
        raise TuneRefused(f"{out_dir} exists and is not empty: a tuning run never overwrites another")
    warnings = list(warnings or [])
    ctxs = [_prepare(it, k, backend, warnings, current) for k, it in enumerate(items)]
    _same_base(ctxs)
    if not any(len(c["objects"]) for c in ctxs):
        raise TuneRefused("no truth point lies in an object of its run: nothing can be scored (points need the "
                          "run's labels under them; docs/TUNING.md)")
    sets = parameter_sets(grid)
    groups: dict[str, list[int]] = {}
    for i, params in enumerate(sets):
        key = group_key(params) if shortcut else str(i)
        groups.setdefault(key, []).append(i)

    scores: list[list[dict | None]] = [[None] * len(ctxs) for _ in sets]
    shared: list[dict | None] = [None] * len(sets)
    n_done = 0
    for k, ctx in enumerate(ctxs):
        labels = _read_labels(ctx["result_dir"], ctx["image_id"])     # one image's labels and pixels at a time
        stack = _read_stack(ctx, labels)
        for members in groups.values():
            if shortcut and len(members) > 1:
                base = {p: v for p, v in sets[members[0]].items() if p != EFFECT_SIZE_KEY}
                base[EFFECT_SIZE_KEY] = 0.0
                t = time.time()
                ps0, pn0 = _detect(ctx, stack, labels, _section(ctx["base"], base, backend))
                dt = time.time() - t
                n_done += 1
                log.info("tune: %s (%d/%d), shared %s run %s: %d spots before the effect-size filter (%.1f s)",
                         ctx["image_id"], k + 1, len(ctxs), backend, json.dumps(base, sort_keys=True), len(ps0), dt)
                for i in members:
                    e = float(sets[i][EFFECT_SIZE_KEY])
                    ps, pn = filter_effect_size(ps0, pn0, e)
                    if on_tables is not None:
                        on_tables(i, k, ps, pn)
                    scores[i][k] = {**_score(ctx, ps), "detect_s": round(dt, 2)}
                    shared[i] = base
            else:
                for i in members:
                    t = time.time()
                    ps, pn = _detect(ctx, stack, labels, _section(ctx["base"], sets[i], backend))
                    dt = time.time() - t
                    n_done += 1
                    if on_tables is not None:
                        on_tables(i, k, ps, pn)
                    scores[i][k] = {**_score(ctx, ps), "detect_s": round(dt, 2)}
                    log.info("tune: %s (%d/%d), set %d of %d %s: %d spots (%.1f s)", ctx["image_id"], k + 1,
                             len(ctxs), i + 1, len(sets), json.dumps(sets[i], sort_keys=True), len(ps), dt)
        del stack, labels

    trials = [_trial(i, sets[i], scores[i], ctxs, backend, split, shared[i], (primary, tie)) for i in range(len(sets))]
    order = sorted(range(len(sets)), key=lambda i: _rank_key(trials[i]), reverse=True)
    for r, i in enumerate(order, start=1):
        trials[i]["rank"] = r
    ranked = [_ranked_row(trials[i], grid) for i in order]
    best_i = order[0]
    best = None
    if any(math.isfinite(v) for v in _rank_key(trials[best_i])):
        best = {"objective": objective_text((primary, tie)), "backend": backend, "keys": dict(sets[best_i])}
    else:
        warnings.append("no parameter set could be scored (fewer than two annotated objects and no matched point): "
                        "no best.yaml was written")
    filtered = shortcut and any(len(m) > 1 for m in groups.values())      # a set was filtered from a shared run
    run = _run_record(run_id, ctxs, grid, backend, split, (primary, tie), max_runs, n_done, warnings, root,
                      time.time() - t_start, filtered)
    _write_run(out_dir, grid, trials, ranked, best, run)
    return {"run_dir": str(out_dir), "best": best, "ranked": ranked, "n_runs": n_done, "warnings": warnings}


# ------------------------------------------------------------------------------ one finished run
def _prepare(item, k: int, backend: str, warnings: list, current) -> dict:
    """What a search needs of one item, read before any pixel of its image: the done marker, the config the
    run was made with, the geometry, the shape of the labels, the channel of the spots stage, the spacing the
    stage used, and the truth points on the run's grid with their objects in the run's labels (the labels
    themselves are read again when the item's detections run, so one image's labels are in memory at a time)."""
    from germquant import provenance
    from germquant.benchmark.spots import label_at_points

    if not isinstance(item, dict):
        raise TuneRefused(f"item {k}: a mapping with image_id, image, result_dir and points, got {item!r}")
    missing = [key for key in ("image_id", "image", "result_dir", "points") if key not in item]
    if missing:
        raise TuneRefused(f"item {k}: missing {missing}")
    image_id = str(item["image_id"])
    image = Path(os.fspath(item["image"]))
    result_dir = Path(os.fspath(item["result_dir"]))
    marker = provenance.read_done_marker(result_dir, image_id)
    if marker is None:
        raise TuneRefused(f"{image_id}: no finished run in {result_dir} (no {image_id}__done.json); run the image "
                          f"first, then annotate it on that run")
    if marker.get("t_index") is not None:
        raise TuneRefused(f"{image_id}: {result_dir} holds one timepoint of a time series; spot tuning takes fixed "
                          f"images in this release")
    manifest, manifest_path = _run_manifest(result_dir)
    if manifest is None or not isinstance(manifest.get("config"), dict):
        raise TuneRefused(f"{image_id}: no run_manifest.json with the run's config in {result_dir} or a folder "
                          f"above it")
    if manifest.get("config_hash") != marker.get("config_hash"):
        raise TuneRefused(f"{image_id}: {manifest_path} records config {manifest.get('config_hash')}, the run's "
                          f"done marker {marker.get('config_hash')}: the manifest is not this run's")
    if "spots" not in (marker.get("enabled_stages") or []) and "spots" not in (manifest.get("enabled_stages") or []):
        warnings.append(f"{image_id}: its run had the spots stage off; tuning runs it on the run's labels")
    base = copy.deepcopy(manifest["config"].get("spots") or {})
    # a single-plane run records itself in its done marker ({"planar": true}, the per-image record); its own
    # run_manifest.json says so too, but `germquant batch` writes one manifest for the whole folder, without
    # per-image geometry
    planar = bool(marker.get("planar") or manifest.get("planar"))
    labels = _read_labels(result_dir, image_id)
    summary = _image_summary(result_dir, image_id)
    ran_on = _summary_file_name(summary)
    if ran_on is not None and ran_on.casefold() != image.name.casefold():
        warnings.append(f"{image_id}: its run read a file named {ran_on}, and the image given is {image.name}; "
                        f"only the image's shape is checked, so make sure it is the same image")
    roles = _roles(summary, image_id)
    if roles.get("foci") is None:
        raise TuneRefused(f"{image_id}: its run bound no channel to the foci role, the channel the spots stage "
                          f"reads")
    spacing = _run_spacing(summary, image_id, planar)
    xy_stride = int(marker.get("xy_stride") or 1)
    z_range = tuple(int(v) for v in marker["z_range"]) if marker.get("z_range") else None
    if planar:
        _check_planar_support(image_id)
        _check_planar_backend(image_id, backend)
    pts = item["points"]
    try:
        if isinstance(pts, pd.DataFrame):
            points_path = None
            pts = check_points(pts.copy(), f"the points of {image_id}")
        else:
            points_path = Path(os.fspath(pts))
            pts = read_points(points_path)
    except FileNotFoundError:
        raise
    except (ValueError, OSError) as e:      # a malformed cell (pandas' ParserError is a ValueError), a locked file
        raise TuneRefused(f"{image_id}: {e}") from None
    frame = to_run_frame(pts, shape=labels.shape, spacing=spacing, xy_stride=xy_stride, z_range=z_range)
    inside = frame["inside"]
    obj = np.zeros(len(pts), np.int64)
    if inside.any():
        obj[inside] = label_at_points(labels, frame["zr"][inside], frame["yr"][inside], frame["xr"][inside])
    changed = int(((obj != pts["object_id"].to_numpy()) & inside).sum())
    if changed:
        warnings.append(f"{image_id}: {changed} truth point(s) lie in another object of results/{image_id} than "
                        f"their file records (the points were saved on other labels); the run's labels are used")
    scored = inside & (obj > 0)
    flags = pts["projected"].to_numpy(dtype=bool)[scored]
    if flags.any() and not flags.all():
        raise TuneRefused(f"{image_id}: its truth points mix points placed on a projection with 3D points; keep "
                          f"one kind per file")
    objects = np.unique(obj[scored])
    if not len(objects):
        warnings.append(f"{image_id}: none of its {len(pts)} truth point(s) lies in an object of its run")
    truth_n = np.array([int((obj[scored] == o).sum()) for o in objects], np.int64)
    truth_xyz = np.column_stack([frame["x_um"][scored], frame["y_um"][scored], frame["z_um"][scored]])
    if current is not None and marker.get("config_hash") != current.get("config_hash"):
        same_spots = json.dumps(base, sort_keys=True) == json.dumps(current.get("spots") or {}, sort_keys=True)
        warnings.append(f"{image_id}: results/{image_id} was made with config {marker.get('config_hash')}, the "
                        f"recipe is now {current.get('config_hash')}"
                        + ("" if same_spots else " with other spots: settings")
                        + "; the search starts from the run's own labels and spot settings")
    return {
        "k": k, "image_id": image_id, "image": image, "result_dir": result_dir, "marker": marker,
        "manifest_path": manifest_path, "config_hash": marker.get("config_hash"), "base": base, "planar": planar,
        "labels_shape": tuple(labels.shape), "roles": roles, "spacing": spacing, "xy_stride": xy_stride, "z_range": z_range,
        "points_path": points_path, "points_sha256": _lf_sha256(points_path) if points_path else None,
        "n_points": int(len(pts)), "n_outside_run": int((~inside).sum()),
        "n_outside_objects": int((inside & (obj == 0)).sum()), "objects": objects, "truth_n": truth_n,
        "truth_xyz": truth_xyz, "projected": bool(flags.any()), "backend": backend,
    }


def _run_manifest(result_dir: Path) -> tuple[dict | None, Path | None]:
    """run_manifest.json of a results folder, or of the nearest folder above it that has one: `germquant batch`
    writes one manifest for its whole output folder, also for the images of subfolders (<out>/<rel>/<id>/).
    The caller checks that the manifest's config hash is the run's."""
    folder = Path(os.path.abspath(result_dir))
    for p in (d / "run_manifest.json" for d in (folder, *folder.parents)):
        try:
            with open(long_path(p), encoding="utf-8") as fh:
                return json.load(fh), p
        except (OSError, ValueError):
            continue
    return None, None


def _read_labels(result_dir: Path, image_id: str) -> np.ndarray:
    import tifffile

    p = result_dir / f"{image_id}__nuclei_labels.tif"
    if not os.path.isfile(long_path(p)):
        raise TuneRefused(f"{image_id}: no label image {p.name} in {result_dir}; tuning needs the run's labels "
                          f"(output.write_label_images must stay on)")
    lab = np.asarray(tifffile.imread(long_path(p)))
    if lab.ndim == 2:
        lab = lab[None]
    if lab.ndim != 3:
        raise TuneRefused(f"{image_id}: {p.name} is not a (Z, Y, X) label image (shape {lab.shape})")
    return lab


def _image_summary(result_dir: Path, image_id: str) -> pd.DataFrame:
    for suffix, reader in ((".csv", pd.read_csv), (".parquet", pd.read_parquet)):
        p = result_dir / f"{image_id}__image_summary{suffix}"
        if os.path.isfile(long_path(p)):
            df = reader(long_path(p))
            if len(df):
                return df
    raise TuneRefused(f"{image_id}: no image summary table in {result_dir}")


def _summary_file_name(summary: pd.DataFrame) -> str | None:
    """The name of the file the run read (the image summary's file_path), or None when it is not recorded. A
    root move keeps the name; another name means another file, or a renamed one."""
    try:
        value = summary["file_path"].iloc[0]
    except (KeyError, IndexError):
        return None
    if not isinstance(value, str) or not value.strip():
        return None
    return re.split(r"[\\/]", value.strip())[-1] or None


def _roles(summary: pd.DataFrame, image_id: str) -> dict:
    try:
        roles = json.loads(str(summary["channel_map_json"].iloc[0]))
    except (KeyError, ValueError, TypeError):
        raise TuneRefused(f"{image_id}: its image summary records no channel roles (channel_map_json)") from None
    return {str(r): (None if i is None else int(i)) for r, i in roles.items()}


def _run_spacing(summary: pd.DataFrame, image_id: str, planar: bool) -> tuple[float, float, float]:
    """The voxel size the run's spots stage used, as its tables record it (an io.voxel_size_um override and
    the stride included); a single-plane (planar) run records no dz and runs its spots with PLANAR_DZ."""
    try:
        dz, dy, dx = (float(summary[c].iloc[0]) for c in ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um"))
    except (KeyError, ValueError, TypeError):
        raise TuneRefused(f"{image_id}: its image summary records no voxel size") from None
    if planar or not math.isfinite(dz):
        from germquant import planar as planar_mod

        dz = float(getattr(planar_mod, "PLANAR_DZ", 1.0))
    if not all(math.isfinite(v) and v > 0 for v in (dz, dy, dx)):
        raise TuneRefused(f"{image_id}: its run records the voxel size {(dz, dy, dx)}, which spot detection cannot use")
    return dz, dy, dx


def _check_planar_support(image_id: str) -> None:
    """A single-plane (planar) run needs the planar spot branch of package P2b (detect_spots(planar=True))."""
    import inspect

    from germquant.spots import detect as D

    if "planar" not in inspect.signature(D.detect_spots).parameters:
        raise TuneRefused(f"{image_id}: its run is single-plane (planar); tuning it needs the planar spot branch "
                          f"of package P2b, which this engine does not have yet")


def _check_planar_backend(image_id: str, backend: str) -> None:
    """A single-plane run detects spots on one plane, which only a backend that declares planar_ok can do (the
    pipeline refuses the others, planar.check_planar); refused here before any detection."""
    from germquant.planar import planar_ok
    from germquant.spots.backends import get_backend

    if not planar_ok(get_backend(backend)):
        raise TuneRefused(f"{image_id}: its run is single-plane (planar), and the {backend} spot backend does not "
                          f"declare planar_ok, so it cannot detect spots on one plane (docs/PLANAR.md)")


def _check_backend(backend: str) -> None:
    """An unknown spot backend is refused before any item is read."""
    from germquant.spots.backends import get_backend

    try:
        get_backend(backend)
    except (TypeError, ValueError) as e:
        raise TuneRefused(str(e)) from None


def _same_base(ctxs: list[dict]) -> None:
    """Every item must come from runs made with one spots: section, so one parameter set means one recipe."""
    first = json.dumps(ctxs[0]["base"], sort_keys=True)
    for c in ctxs[1:]:
        if json.dumps(c["base"], sort_keys=True) != first:
            raise TuneRefused(f"{ctxs[0]['image_id']} and {c['image_id']} were run with different spots: settings; "
                              f"run them with one recipe before tuning them together")


def _read_stack(ctx: dict, labels: np.ndarray):
    """The image at the run's geometry, through the reader process_image uses (pipeline._read_input), checked
    against the run's labels."""
    from germquant import pipeline

    _, stack = pipeline._read_input(ctx["image"], xy_stride=ctx["xy_stride"], z_range=ctx["z_range"])
    foci = ctx["roles"]["foci"]
    if foci >= stack.n_channels:
        raise ValueError(f"{ctx['image_id']}: {ctx['image']} has {stack.n_channels} channel(s); its run read the "
                         f"foci channel at index {foci}")
    shape = tuple(stack.data[foci].shape)
    if shape != tuple(labels.shape) or shape != ctx["labels_shape"]:
        raise ValueError(f"{ctx['image_id']}: {ctx['image']} reads as {shape} at the run's geometry, its labels are "
                         f"{tuple(labels.shape)}: the image is not the one the run read")
    return stack


def _section(base: dict, params: dict, backend: str) -> dict:
    """The spots: section of a run with ``params`` ({dotted spots key: value}) laid over ``base``."""
    section = copy.deepcopy(base)
    section["backend"] = backend
    for key, value in params.items():
        parts = key.split(".")[1:]
        node = section
        for part in parts[:-1]:
            nxt = node.get(part)
            if not isinstance(nxt, dict):
                nxt = {}
                node[part] = nxt
            node = nxt
        node[parts[-1]] = value
    return section


def stage_kwargs(section: dict) -> dict:
    """The keyword arguments the spots stage of process_image passes to detect_spots for a config whose
    spots: section is ``section`` (pipeline.process_image, _spots): every SpotMAX parameter with the stage's
    own code defaults, the backend and the raw section. tests/test_tune_spots.py checks that a detection with
    these arguments reproduces a run's spots table."""
    g = section.get
    return dict(
        spot_radius_um=float(g("spot_radius_um", 0.3)),
        gauss_sigma_um=float(g("gauss_sigma_um", 0.08)),
        thresholding_method=g("thresholding_method", "threshold_triangle"),
        effect_size_metric=g("effect_size_metric", "spot_vs_backgr_effect_size_glass"),
        effect_size_min=float(g("effect_size_min", 3.0)),
        merge_z_columns=bool(g("merge_z_columns", True)),
        z_merge_gap_um=float(g("z_merge_gap_um", 0.8)),
        z_merge_valley_frac=float(g("z_merge_valley_frac", 0.8)),
        max_spot_candidates=int(g("max_spot_candidates", 30000)),
        backend=g("backend", "spotmax"),
        cfg_section=dict(section),
    )


def _detect(ctx: dict, stack, labels: np.ndarray, section: dict):
    """One detection on the run's labels with the spots: section ``section``, as the spots stage makes it, but
    for the tables' marker column: the stage writes the channel map's marker name there, which a finished run
    keeps only in the rows of its spots table, so the tuner passes the role name foci. The marker only labels
    the tables (the backend contract of germquant.spots.backends), and scoring does not read it."""
    from germquant import pipeline
    from germquant.spots import detect as D

    kw = stage_kwargs(section)
    extra = pipeline.spot_backend_channels(kw["backend"], ctx["roles"], stack)   # {} for spotmax and log
    if ctx["planar"]:
        extra["planar"] = True
    return D.detect_spots(stack.data[ctx["roles"]["foci"]], labels, ctx["spacing"], marker="foci", **kw, **extra)


# ------------------------------------------------------------------------------ scores
def _score(ctx: dict, per_spot: pd.DataFrame) -> dict:
    """One detection scored inside the annotated objects of one image."""
    from germquant.benchmark.spots import match_spots
    from germquant.benchmark.stats import precision_recall_f1

    objects = ctx["objects"]
    nuc = per_spot["nucleus_id"].to_numpy(dtype=np.int64) if len(per_spot) else np.zeros(0, np.int64)
    ours = per_spot.loc[np.isin(nuc, objects)] if len(per_spot) else per_spot
    if len(ours):
        counts = pd.Series(ours["nucleus_id"].to_numpy(dtype=np.int64)).value_counts()
        ours_n = np.array([int(counts.get(int(o), 0)) for o in objects], np.int64)
        ours_xyz = ours[["x_um", "y_um", "z_um"]].to_numpy(dtype=float)
    else:
        ours_n = np.zeros(len(objects), np.int64)
        ours_xyz = np.zeros((0, 3))
    truth_xyz = ctx["truth_xyz"]
    if ctx["projected"]:                       # points placed on a projection are matched in 2D (design 6.7)
        ours_xyz = ours_xyz.copy()
        truth_xyz = truth_xyz.copy()
        ours_xyz[:, 2] = 0.0
        truth_xyz[:, 2] = 0.0
    matched = match_spots(ours_xyz, truth_xyz, tol_xy_um=TOL_XY_UM, tol_z_um=TOL_Z_UM).n_matched
    det = precision_recall_f1(matched, len(ours_xyz), len(truth_xyz))
    return {"truth_n": ctx["truth_n"], "ours_n": ours_n, "n_matched": int(matched), "n_ours": int(len(ours_xyz)),
            "n_truth": int(len(truth_xyz)), "object_ccc": _ccc(ctx["truth_n"], ours_n),
            "point_precision": det["precision"], "point_recall": det["recall"], "point_f1": det["f1"],
            "n_spots_all": int(len(per_spot))}


def _ccc(truth, ours) -> float:
    """Lin's CCC of per-object counts (germquant.benchmark.stats.lin_ccc), and 1.0 for two equal series of at
    least two objects: lin_ccc is NaN there when every object holds the same count (0 / 0), and the set that
    finds exactly the annotated counts must not rank below every other set. With one object it stays NaN (one
    pair defines no concordance), and against a constant truth every other series scores 0.0."""
    from germquant.benchmark.stats import lin_ccc

    t, o = np.asarray(truth, float), np.asarray(ours, float)
    if len(t) >= 2 and t.shape == o.shape and np.array_equal(t, o):
        return 1.0
    return lin_ccc(t, o)


def _num(v):
    """A finite float, or None (JSON has no NaN)."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _trial(i: int, params: dict, per_item: list[dict], ctxs: list[dict], backend: str, split: str, shared,
           objective: tuple[str, str]) -> dict:
    from germquant.benchmark.stats import precision_recall_f1

    truth = np.concatenate([s["truth_n"] for s in per_item]) if per_item else np.zeros(0)
    ours = np.concatenate([s["ours_n"] for s in per_item]) if per_item else np.zeros(0)
    m = sum(s["n_matched"] for s in per_item)
    no = sum(s["n_ours"] for s in per_item)
    nt = sum(s["n_truth"] for s in per_item)
    det = precision_recall_f1(m, no, nt)
    pooled = {"object_ccc": _num(_ccc(truth, ours)), "point_f1": _num(det["f1"]),
              "point_precision": _num(det["precision"]), "point_recall": _num(det["recall"]),
              "n_objects": int(len(truth)), "n_points": int(nt), "n_matched": int(m), "n_ours": int(no),
              "truth_count": int(truth.sum()), "ours_count": int(ours.sum())}
    images = []
    for s, c in zip(per_item, ctxs):
        images.append({"image_id": c["image_id"], "object_ccc": _num(s["object_ccc"]), "point_f1": _num(s["point_f1"]),
                       "point_precision": _num(s["point_precision"]), "point_recall": _num(s["point_recall"]),
                       "n_objects": int(len(s["truth_n"])), "n_points": s["n_truth"], "n_matched": s["n_matched"],
                       "n_ours": s["n_ours"], "n_spots_all": s["n_spots_all"],
                       "n_points_outside_objects": c["n_outside_objects"], "n_points_outside_run": c["n_outside_run"],
                       "detect_s": s["detect_s"]})
    primary, tie = objective
    pv, tv = pooled[primary], pooled[tie]
    return {"trial": i, "params": dict(params), "backend": backend, "split": split, "shared_detection": shared,
            "objective": {"name": primary, "value": None if pv is None else round(pv, CCC_DECIMALS),
                          "tie_break": tie, "tie_break_value": tv},
            "pooled": pooled, "images": images}


def _rank_key(trial: dict) -> tuple[float, float]:
    o = trial["objective"]
    return (-math.inf if o["value"] is None else float(o["value"]),
            -math.inf if o["tie_break_value"] is None else float(o["tie_break_value"]))


def _ranked_row(trial: dict, grid: dict) -> dict:
    o = trial["objective"]
    row = {"rank": trial["rank"], "objective": o["name"], "objective_value": o["value"], "tie_break": o["tie_break"],
           "tie_break_value": None if o["tie_break_value"] is None else round(o["tie_break_value"], CCC_DECIMALS)}
    for key in grid:
        row[key] = trial["params"][key]
    return row


# ------------------------------------------------------------------------------ the run folder
def _lf_sha256(path) -> str | None:
    """sha256 of a text file with CRLF turned into LF (as germquant.packs.format.lf_sha256), so a points file
    re-saved by Excel and a done marker written on Windows give the digest they give on Linux."""
    import hashlib

    if path is None:
        return None
    with open(long_path(path), "rb") as fh:
        return hashlib.sha256(fh.read().replace(b"\r\n", b"\n")).hexdigest()


def _rel(path, root) -> str:
    """A path for run.json: relative to the project with forward slashes when it lies inside it."""
    if path is None:
        return None
    p = Path(path)
    if root is not None:
        try:
            return p.resolve().relative_to(Path(root).resolve()).as_posix()
        except ValueError:
            pass
    return str(p)


def _run_record(run_id, ctxs, grid, backend, split, objective, max_runs, n_done, warnings, root, elapsed,
                filtered) -> dict:
    from germquant import __version__, provenance

    items = []
    for c in ctxs:
        items.append({"image_id": c["image_id"], "image": _rel(c["image"], root), "result_dir": _rel(c["result_dir"], root),
                      "config_hash": c["config_hash"],
                      "done_sha256": _lf_sha256(c["result_dir"] / f"{c['image_id']}__done.json"),
                      "xy_stride": c["xy_stride"], "z_range": list(c["z_range"]) if c["z_range"] else None,
                      "spacing_um": list(c["spacing"]), "planar": c["planar"],
                      "points": _rel(c["points_path"], root), "points_sha256": c["points_sha256"],
                      "n_points": c["n_points"], "n_points_scored": int(len(c["truth_xyz"])),
                      "n_points_outside_objects": c["n_outside_objects"], "n_points_outside_run": c["n_outside_run"],
                      "n_objects": int(len(c["objects"])), "projected": c["projected"]})
    return {"format": 1, "run_id": run_id, "created": provenance.run_timestamp(), "germquant": __version__,
            "git_sha": provenance.git_sha(), "backend": backend, "split": split,
            "objective": objective_text(objective), "objective_metrics": list(objective),
            "objective_decimals": CCC_DECIMALS,
            "shortcut": "run once, then filter" if filtered else "one detection per set",
            "match": {"tol_xy_um": TOL_XY_UM, "tol_z_um": TOL_Z_UM, "scored": "annotated objects"},
            "grid": grid, "n_sets": len(parameter_sets(grid)), "n_runs": n_done, "max_runs": max_runs,
            "base_spots": ctxs[0]["base"], "items": items, "warnings": list(warnings), "elapsed_s": round(elapsed, 1)}


def _write_run(out_dir: Path, grid: dict, trials: list[dict], ranked: list[dict], best: dict | None, run: dict) -> None:
    import yaml

    os.makedirs(long_path(out_dir), exist_ok=True)
    write_text_atomic(out_dir / FILES["grid"], yaml.safe_dump(grid, sort_keys=False, allow_unicode=False,
                                                               default_flow_style=None))
    lines = [json.dumps(t, sort_keys=False, ensure_ascii=True, allow_nan=False) for t in trials]
    write_text_atomic(out_dir / FILES["trials"], "".join(line + "\n" for line in lines))
    cols = ["rank", "objective", "objective_value", "tie_break", "tie_break_value", *grid]
    write_text_atomic(out_dir / FILES["ranked"], pd.DataFrame(ranked, columns=cols).to_csv(index=False,
                                                                                             lineterminator="\n"))
    if best is not None:
        write_text_atomic(out_dir / FILES["best"], yaml.safe_dump(best, sort_keys=False, allow_unicode=False))
    write_text_atomic(out_dir / FILES["run"], json.dumps(run, indent=1, ensure_ascii=True, allow_nan=False) + "\n")


def read_best(run_dir) -> dict:
    """best.yaml of a tuning run, checked: {objective (text), backend (a name), keys {dotted spots key:
    value}}. TuneRefused for a run without one (nothing could be scored) or a malformed file."""
    import yaml

    p = Path(os.fspath(run_dir)) / FILES["best"]
    if not os.path.isfile(long_path(p)):
        raise TuneRefused(f"{p.parent}: no best.yaml (not a tuning run, or no parameter set could be scored)")
    with open(long_path(p), encoding="utf-8-sig") as fh:
        best = yaml.safe_load(fh)
    if not isinstance(best, dict) or set(best) != {"objective", "backend", "keys"} or not isinstance(best["keys"], dict) \
            or not best["keys"] or not isinstance(best["backend"], str):
        raise TuneRefused(f"{p}: best.yaml holds objective, backend and keys {{dotted spots key: value}}")
    check_grid({k: [v] for k, v in best["keys"].items()})
    return best
