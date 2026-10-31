"""Spot preset tuning from truth points (design 6.8 and 10 step 9; plan D24 and D28; docs/TUNING.md).

Spots are tuned, not trained: a search runs the spot detection again on the labels of finished runs, once per
parameter set of a grid of spot keys, scores each set against truth points and ranks the sets by a fixed
objective (the per-object Lin's CCC over the annotated objects, rounded to 3 decimals, then the point F1).

* ``germquant.tune.points``: the truth-point files, annotations/<image_id>__points_<annotator>.csv
  (``write_points``, ``read_points``);
* ``germquant.tune.grid``: grids (``default_grid``, ``load_grid``, ``n_runs``);
* ``germquant.tune.effect_filter``: SpotMAX's effect-size filter (``filter_effect_size``);
* ``germquant.tune.spots``: the tuners (``tune_spotmax``, ``tune_grid``) and the project entry
  ``tune_spots``;
* ``germquant.tune.apply``: ``apply_tuning``, which writes the best set into a project's recipe.yaml.

Importing germquant.tune loads no torch, CuPy, SpotMAX, Cellpose, bioio or dask.array (plan D28): a tuner
imports germquant.pipeline (torch first) when it runs, before any reader or backend. Of the benchmark and
validation code, tune imports only germquant.benchmark.stats, matching, nuclei3d and spots and
germquant.validate.compare (the seal of design 6.7, enforced by tests/test_tune_seal.py).
"""
from .apply import apply_tuning
from .effect_filter import filter_effect_size
from .grid import GridError, TuneRefused, TuningWarning, default_grid, load_grid, n_runs
from .points import read_points, write_points
from .spots import OBJECTIVE, OBJECTIVE_TEXT, tune_grid, tune_spotmax, tune_spots

__all__ = [
    "OBJECTIVE",
    "OBJECTIVE_TEXT",
    "GridError",
    "TuneRefused",
    "TuningWarning",
    "apply_tuning",
    "default_grid",
    "filter_effect_size",
    "load_grid",
    "n_runs",
    "read_points",
    "tune_grid",
    "tune_spotmax",
    "tune_spots",
    "write_points",
]
