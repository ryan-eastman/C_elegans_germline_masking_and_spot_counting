"""Collecting a project's results (docs/PROJECTS.md, "Collecting the tables"; design 7.8 and 3.12; P5 acceptance 14).

``collect_project(project)`` stacks the per-image tables under results/ into results/batch_<table>.csv and rebuilds
results/batch_summary.csv with batch.collect, the code `germquant collect` and `germquant batch` use, for the
project's own images only (its image ids; an image with timepoints by its timepoint ids <id>_t<NNN>), and with the
recipe's qc.germline_outlier_factor. The plain `germquant collect` has no recipe and keeps the factor 1.8, which
would flag dense fields of the outside packs as framing outliers (their base sets 1.0e+9 to switch it off).
Returns ``{"root", "results", "outlier_factor", "images", "counts", "summary"}``; ``counts`` is batch.collect's
``{table: rows}``.
"""
from __future__ import annotations

import os
from pathlib import Path

from ..fsutil import long_path
from .manifest import ProjectError


def collect_ids(manifest: dict) -> set[str]:
    """The ids collect stacks: every image of the project not excluded, a series by its timepoint ids."""
    from .run import is_series

    ids: set[str] = set()
    for e in manifest.get("images") or []:
        if e.get("excluded"):
            continue
        iid = str(e["id"]).strip()
        if is_series(e):
            ids.update(f"{iid}_t{k:03d}" for k in range(int(e["timepoints"])))
        else:
            ids.add(iid)
    return ids


def collect_project(project) -> dict:
    """Collect the per-image tables of results/ with the recipe's qc.germline_outlier_factor (module docstring).
    ProjectError when the project has no results yet."""
    from .. import batch
    from ..packs import use_project
    from . import Project, _root, project_config
    from .manifest import read_manifest
    from .run import RESULTS_DIR

    root = _root(project)
    m = read_manifest(root)
    with use_project(root):
        cfg = project_config(Project(root=root, manifest=m))
    factor = float(cfg.get("qc.germline_outlier_factor", 1.8))
    results = Path(root) / RESULTS_DIR
    if not os.path.isdir(long_path(results)):
        raise ProjectError(f"{root} has no results yet ({results} is missing): run the project first with "
                           f"`germquant project run DIR`")
    ids = collect_ids(m)
    counts = batch.collect(results, image_ids=ids, outlier_factor=factor)
    return {"root": str(root), "results": str(results), "outlier_factor": factor, "images": len(ids),
            "counts": counts, "summary": str(results / "batch_summary.csv")}
