"""Tuning grids: which spot keys a tuning run searches, and over which values (design 6.8; D24).

A grid is a mapping of dotted ``spots.`` keys to the list of values to try, for example::

    spots.log.min_snr: [1.5, 2, 2.5, 3, 4, 5]
    spots.log.diameter_um: [0.4, 0.5]

Every combination is one parameter set (the product, in the order of the keys and of their values); every
other key keeps the value of the recipe the runs were made with. ``load_grid`` reads such a YAML file (the
``--grid`` file of ``germquant tune-spots``, and the grid.yaml of every tuning run, so a run can be
repeated); ``check_grid`` validates one in memory.

The keys a known backend reads are checked by name and type: for ``spotmax`` the SpotMAX parameters of the
spots stage (``spots.detect.SPOTMAX_KEYS``) and the optional keys of a pack recipe
(``SPOTMAX_OPTIONAL_KEYS``); for ``log`` the ``spots.log.*`` keys and the shared z-merge and flood-guard
keys. A number is stored with the type of its code default (``2`` becomes ``2.0`` for a float key), so
best.yaml and recipe.yaml hold the same values a run reads. A backend from an entry point takes any
``spots.`` key with a scalar value. ``spots.backend``, ``spots.enabled`` and ``spots.instances`` are never
tuned.

Default grids, used when the recipe declares none (no pack ships a tuning grid yet; engine_specs 10.14):
``spots.effect_size_min`` in {0, 0.5, 1, 1.5, 2, 3, 4, 5} for spotmax and ``spots.log.min_snr`` in
{1.5, 2, 2.5, 3, 4, 5} for log. Any other backend needs a grid file.

``n_runs`` counts the detections a search makes: SpotMAX runs once per group of parameter sets that differ
only in ``spots.effect_size_min`` (run once, then filtered; design 6.8), every other backend once per
parameter set, each on every image.
"""
from __future__ import annotations

import itertools
import json
import math
import os
from pathlib import Path

EFFECT_SIZE_KEY = "spots.effect_size_min"
RESERVED_KEYS = ("spots.backend", "spots.enabled", "spots.instances")
DEFAULT_GRIDS = {
    "spotmax": {EFFECT_SIZE_KEY: [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]},
    "log": {"spots.log.min_snr": [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]},
}


class TuneRefused(ValueError):
    """A tuning run refused before any detection: a grid that cannot be searched, too many runs for
    --max-runs, an image without a finished run or without truth points (``germquant tune-spots`` exits 2)."""


class GridError(TuneRefused):
    """A grid that cannot be searched: an unknown or reserved key, an empty list, a value of the wrong type."""


class TuningWarning(UserWarning):
    """Something to check before a tuning result is used, such as a recipe that changed after the tuning run
    (apply_tuning warns with it; ``germquant tune-spots --apply`` prints it as a WARNING line)."""


def _yaml_problem(e: BaseException) -> str:
    """A YAML parser's error, or the UnicodeDecodeError of a file that is not UTF-8 text, on one line: the
    problem and where it is, so that a refusal is one ERROR line."""
    mark = getattr(e, "problem_mark", None)
    if getattr(e, "problem", None) and mark is not None:
        return f"{e.problem} at line {mark.line + 1}, column {mark.column + 1}"
    return " ".join(str(e).split())


def _key_types(backend: str) -> dict[str, type | tuple] | None:
    """The tunable keys of a built-in backend with the type of their code default, or None for a backend
    whose keys are not known here (an entry point)."""
    if backend == "spotmax":
        from germquant.spots.detect import SPOTMAX_DEFAULTS, SPOTMAX_OPTIONAL_KEYS

        types = {f"spots.{k}": type(v) for k, v in SPOTMAX_DEFAULTS.items()}
        types.update({"spots.spot_radius_z_um": float, "spots.remove_hot_pixels": bool, "spots.sharpen": bool})
        missing = [k for k in SPOTMAX_OPTIONAL_KEYS if f"spots.{k}" not in types]
        if missing:                                   # a new optional key must be typed here before it is tuned
            raise RuntimeError(f"germquant.tune.grid has no type for the optional SpotMAX keys {missing}")
        return types
    if backend == "log":
        from germquant.spots.log_detect import LOG_DEFAULTS, SHARED_DEFAULTS

        types = {f"spots.log.{k}": float for k in LOG_DEFAULTS}       # the None defaults are derived lengths
        types.update({f"spots.{k}": type(v) for k, v in SHARED_DEFAULTS.items()})
        return types
    return None


def _coerce(key: str, value, want):
    """``value`` as the type ``want`` of the key's code default, or GridError naming the key."""
    if isinstance(value, (dict, list, tuple, set)) or value is None:
        raise GridError(f"{key}: every value must be a scalar, got {value!r}")
    if isinstance(value, str) and not value.isascii():
        raise GridError(f"{key}: values must be ASCII text, got {value!r}")
    if want is None:
        if isinstance(value, float) and not math.isfinite(value):
            raise GridError(f"{key}: values must be finite numbers, got {value!r}")
        return value
    if want is bool:
        if not isinstance(value, bool):
            raise GridError(f"{key} takes true or false, got {value!r}")
        return value
    if want is int:
        if isinstance(value, bool) or not isinstance(value, int):
            raise GridError(f"{key} takes whole numbers, got {value!r}")
        return int(value)
    if want is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise GridError(f"{key} takes finite numbers, got {value!r}")
        return float(value)
    if want is str:
        if not isinstance(value, str):
            raise GridError(f"{key} takes text, got {value!r}")
        return value
    raise GridError(f"{key}: no rule for values of type {want.__name__}")   # pragma: no cover


def check_grid(grid, backend: str | None = None) -> dict:
    """A validated copy of ``grid`` ({dotted spots key: [values]}, keys in their given order). With
    ``backend`` a built-in (spotmax, log) the keys must be ones that backend reads, and numbers take the type
    of the key's code default. GridError for an empty grid, a key outside ``spots.``, a reserved key, an
    empty or repeated list of values, or a value of the wrong type."""
    if not isinstance(grid, dict) or not grid:
        raise GridError("a grid is a non-empty mapping of dotted spots keys to lists of values")
    types = _key_types(backend) if backend else None
    out: dict = {}
    for key, values in grid.items():
        if not isinstance(key, str) or not key.startswith("spots.") or not key.isascii() or key.endswith("."):
            raise GridError(f"grid key {key!r}: tuning searches spot keys only, written as spots.<key>")
        if key in RESERVED_KEYS or any(key.startswith(r + ".") for r in RESERVED_KEYS):
            raise GridError(f"grid key {key!r} is never tuned (the backend, the switch and the instances are "
                            f"recipe choices)")
        if types is not None and key not in types:
            raise GridError(f"grid key {key!r} is not a key the {backend} backend reads; it reads: "
                            f"{', '.join(sorted(types))}")
        if not isinstance(values, (list, tuple)) or not values:
            raise GridError(f"{key}: give a non-empty list of values to try, got {values!r}")
        want = types.get(key) if types is not None else None
        vals = [_coerce(key, v, want) for v in values]
        seen = [json.dumps(v) for v in vals]
        if len(set(seen)) != len(seen):
            raise GridError(f"{key}: a value appears twice in {values!r}")
        out[key] = vals
    return out


def default_grid(cfg) -> dict:
    """The grid tune-spots searches when the recipe declares none: ``spots.effect_size_min`` in {0, 0.5,
    1, 1.5, 2, 3, 4, 5} for spotmax and ``spots.log.min_snr`` in {1.5, 2, 2.5, 3, 4, 5} for log (D24); every
    other key keeps the recipe's value. GridError for another backend, which needs a grid file."""
    backend = str(cfg.get("spots.backend", "spotmax") if hasattr(cfg, "get") else cfg)
    if backend not in DEFAULT_GRIDS:
        raise GridError(f"the {backend} spot backend has no default tuning grid; pass one with --grid FILE "
                        f"(docs/TUNING.md)")
    return {k: list(v) for k, v in DEFAULT_GRIDS[backend].items()}


def load_grid(path) -> dict:
    """The grid in the YAML file ``path`` ({dotted spots key: [values]}), validated by ``check_grid``
    (without a backend: the tuner checks the keys against its backend). GridError for a file that is not
    YAML in UTF-8 (or ASCII) text, or not such a mapping; FileNotFoundError for a missing file."""
    import yaml

    from germquant.fsutil import long_path

    p = Path(os.fspath(path))
    with open(long_path(p), encoding="utf-8-sig") as fh:
        try:
            data = yaml.safe_load(fh)
        except (yaml.YAMLError, UnicodeDecodeError) as e:     # UnicodeDecodeError: not UTF-8 text
            raise GridError(f"{p}: not a YAML file ({_yaml_problem(e)})") from e
    if not isinstance(data, dict):
        raise GridError(f"{p}: a grid file is a mapping of dotted spots keys to lists of values")
    try:
        return check_grid(data)
    except GridError as e:
        raise GridError(f"{p}: {e}") from None


def parameter_sets(grid: dict) -> list[dict]:
    """Every combination of the grid's values, as {key: value}: the product in the order of the keys and of
    their values (the first key varies slowest)."""
    keys = list(grid)
    return [dict(zip(keys, combo)) for combo in itertools.product(*(grid[k] for k in keys))]


def group_key(params: dict) -> str:
    """The identity of a SpotMAX detection apart from its effect-size filter: the parameter set without
    ``spots.effect_size_min``, as canonical JSON. Sets with the same key share one SpotMAX run; sharpen is
    part of the key like every other parameter, so no threshold is ever carried across a sharpen value."""
    return json.dumps({k: v for k, v in params.items() if k != EFFECT_SIZE_KEY}, sort_keys=True)


def n_groups(grid: dict, backend: str) -> int:
    """The number of detections a search makes on one image: the SpotMAX groups (sets that differ only in
    ``spots.effect_size_min`` share one run) for spotmax, else one per parameter set."""
    sizes = [len(v) for k, v in grid.items() if not (backend == "spotmax" and k == EFFECT_SIZE_KEY)]
    return math.prod(sizes) if sizes else 1


def n_runs(grid: dict, n_items: int, backend: str) -> int:
    """The detections a search on ``n_items`` images makes (what --max-runs bounds): ``n_groups`` times the
    number of images."""
    return n_groups(grid, backend) * max(0, int(n_items))
