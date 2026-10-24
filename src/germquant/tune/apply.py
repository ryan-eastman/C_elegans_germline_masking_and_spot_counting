"""Apply a tuning run to a project: write the best parameter set into recipe.yaml (design 10 step 9 and 8.7;
plan D19 and section 9; ``germquant tune-spots --project DIR --apply tuning/<run>``).

The keys of best.yaml go to germquant.project.write_recipe with origin "tuning", so recipe.yaml gains them
on top of the pack's recipe (``extends: pack:<id>@<version>/<recipe>``), the config hash changes and the
project reads validated: false, as after any edit. Only the keys of the run's grid are written (D19: a tuning
edit sets the grid's keys, which may be keys the pack does not declare as tunables, such as
spots.spot_radius_z_um or spots.sharpen), and the run must have tuned the backend the recipe uses now. When
the recipe is no longer the config the run's images were made with (run.json), apply_tuning warns with a
TuningWarning and still writes. Nothing heavy is imported.
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

from germquant.fsutil import long_path

from .grid import TuneRefused, TuningWarning, load_grid
from .spots import FILES, read_best


def apply_tuning(project, run_dir, *, plugin_version=None) -> str:
    """Write the best parameter set of the tuning run ``run_dir`` (relative to the project, for example
    tuning/<run_id>) into the project's recipe.yaml, through write_recipe with origin "tuning", and return the
    project's new config hash. ``project`` is a germquant.project.Project or its folder; ``plugin_version``
    goes into the manifest's saved stamp. TuneRefused for a run without best.yaml, keys outside the run's
    grid, or a run that tuned another backend than the recipe's; a TuningWarning when the recipe changed after
    the run's images were run."""
    open_project, project_config, write_recipe = _project_api()
    if isinstance(project, (str, os.PathLike)):
        project = open_project(project)
    run = Path(os.fspath(run_dir))
    if not run.is_absolute():
        run = Path(project.root) / run
    best = read_best(run)
    grid_file = run / FILES["grid"]
    if not os.path.isfile(long_path(grid_file)):
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
    for message in _changed_since(run, cfg):
        warnings.warn(message, TuningWarning, stacklevel=2)
    return write_recipe(project, edits=dict(best["keys"]), origin="tuning", plugin_version=plugin_version)


def _changed_since(run: Path, cfg) -> list[str]:
    """A warning when the recipe's config is not the one the tuning run's images were run with (the config
    hashes run.json records): the search started from those runs' spot settings, and the best keys are written
    on top of the recipe as it is now. Nothing when run.json is missing or records no config hash."""
    try:
        with open(long_path(run / FILES["run"]), encoding="utf-8") as fh:
            record = json.load(fh)
        made = sorted({str(it["config_hash"]) for it in record["items"] if it.get("config_hash")})
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return []
    now = getattr(cfg, "hash", None)
    other = [h for h in made if h != now]
    if now is None or not other:
        return []
    return [f"{run.name}: its images were run with config {', '.join(other)}, and the recipe is now {now}; the "
            f"best keys were searched on those runs' spot settings and are written on top of the recipe as it "
            f"is now: check its other spot keys, or run the images again and tune again"]


def _project_api():
    """open_project, project_config and write_recipe of germquant.project (a seam the tests replace)."""
    from germquant import project as P

    return P.open_project, P.project_config, P.write_recipe
