"""Apply a tuning run to a project: write the best parameter set into recipe.yaml (design 10 step 9 and 8.7;
plan D19 and section 9; ``germquant tune-spots --project DIR --apply tuning/<run>``).

The keys of best.yaml go to germquant.project.write_recipe with origin "tuning", so recipe.yaml gains them
on top of the pack's recipe (``extends: pack:<id>@<version>/<recipe>``), the config hash changes and the
project reads validated: false, as after any edit. Only the keys of the run's grid are written (D19: a tuning
edit sets the grid's keys, which may be keys the pack does not declare as tunables, such as
spots.spot_radius_z_um or spots.sharpen), and the run must have tuned the backend the recipe uses now.
Nothing heavy is imported.
"""
from __future__ import annotations

import os
from pathlib import Path

from .grid import TuneRefused, load_grid
from .spots import FILES, read_best


def apply_tuning(project, run_dir, *, plugin_version=None) -> str:
    """Write the best parameter set of the tuning run ``run_dir`` (relative to the project, for example
    tuning/<run_id>) into the project's recipe.yaml, through write_recipe with origin "tuning", and return the
    project's new config hash. ``project`` is a germquant.project.Project or its folder; ``plugin_version``
    goes into the manifest's saved stamp. TuneRefused for a run without best.yaml, keys outside the run's
    grid, or a run that tuned another backend than the recipe's."""
    open_project, project_config, write_recipe = _project_api()
    if isinstance(project, (str, os.PathLike)):
        project = open_project(project)
    run = Path(os.fspath(run_dir))
    if not run.is_absolute():
        run = Path(project.root) / run
    best = read_best(run)
    grid_file = run / FILES["grid"]
    if not grid_file.is_file():
        raise TuneRefused(f"{run}: no grid.yaml; apply only takes a folder a tuning run wrote")
    grid = load_grid(grid_file)
    extra = sorted(set(best["keys"]) - set(grid))
    if extra:
        raise TuneRefused(f"{run}: best.yaml sets {extra}, which the run's grid does not search; a tuning edit "
                          f"writes the grid's keys only")
    cfg = project_config(project)
    backend = str(cfg.get("spots.backend", "spotmax"))
    if best["backend"] != backend:
        raise TuneRefused(f"{run} tuned the {best['backend']} backend; the recipe now uses {backend}")
    return write_recipe(project, edits=dict(best["keys"]), origin="tuning", plugin_version=plugin_version)


def _project_api():
    """open_project, project_config and write_recipe of germquant.project (a seam the tests replace)."""
    from germquant import project as P

    return P.open_project, P.project_config, P.write_recipe
