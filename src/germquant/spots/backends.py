"""Pluggable spot detection backends: a registry mapping a backend name to a detector callable.

Contract. A backend is a callable

    backend(image, labels, spacing, *, cfg_section, marker) -> (per_spot, per_nucleus)

where ``image`` is the (Z, Y, X) foci channel, ``labels`` the (Z, Y, X) nucleus label image (0 =
background), ``spacing`` the (dz, dy, dx) voxel size in microns, ``cfg_section`` the raw ``spots:``
mapping of the config with the resolved SpotMAX parameters overlaid (``detect_spots`` builds it; a
backend reads its own keys from it with code defaults and ignores the rest) and ``marker`` the marker
name for the tables. It returns two DataFrames with EXACTLY the columns ``PER_SPOT_COLS`` and
``PER_NUC_COLS`` below (positions in microns, one per-nucleus row for every label, a zero count being a
real zero) and the ``detector`` column of both set to the backend's registry name. ``detect_spots``
checks the columns and warns when the detector column does not carry the name.

Built-ins: ``spotmax`` (the validated default, Imaris cross-validated; its outputs are byte-identical to
the pre-registry code) and ``log`` (an EXPERIMENTAL classical Laplacian-of-Gaussian baseline that is
never selected unless ``spots.backend: log`` is set). Third-party packages add backends through the
entry-point group ``germquant.spot_backends`` (``name = "package.module:function"``); a failing entry
point is logged and skipped, never fatal, and an entry point can not replace a registered name (the
validated default can not be swapped out from outside). docs/BACKENDS.md has the user-facing text.
"""
from __future__ import annotations

import logging
import re
from collections.abc import Callable, Iterable
from typing import Any

import pandas as pd

log = logging.getLogger(__name__)

ENTRY_POINT_GROUP = "germquant.spot_backends"
DEFAULT_BACKEND = "spotmax"

# the table contract (schema.SPOTS for the per-spot table; the per-nucleus table is internal to the stage)
PER_SPOT_COLS = [
    "spot_id", "nucleus_id", "marker", "z_um", "y_um", "x_um",
    "effect_size", "spot_mean_intensity", "detector",
]
PER_NUC_COLS = ["nucleus_id", "marker", "n_spots", "mean_spot_intensity", "detector"]

SpotBackend = Callable[..., tuple[pd.DataFrame, pd.DataFrame]]

_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")   # the name is written into QC flag text
_REGISTRY: dict[str, SpotBackend] = {}
_discovered = False


def register_backend(name: str, fn: SpotBackend | None = None, *, replace: bool = False):
    """Register ``fn`` as the spot backend ``name``. Usable as a decorator (``@register_backend("x")``).
    Re-registering a name raises unless ``replace=True`` (so a typo can not silently shadow a backend)."""
    _check_name(name)

    def _register(f: SpotBackend) -> SpotBackend:
        if not callable(f):
            raise TypeError(f"spot backend {name!r} must be callable, got {type(f).__name__}")
        if name in _REGISTRY and not replace:
            raise ValueError(f"spot backend {name!r} is already registered; pass replace=True to override it")
        _REGISTRY[name] = f
        return f

    return _register if fn is None else _register(fn)


def get_backend(name: str) -> SpotBackend:
    """The backend registered as ``name`` (built-ins, in-process registrations and entry points).
    Raises ValueError naming the known backends for an unknown name."""
    if not isinstance(name, str):
        raise TypeError(f"spot backend name must be a string, got {type(name).__name__}: {name!r}")
    discover_backends()
    try:
        return _REGISTRY[name]
    except KeyError:
        raise ValueError(
            f"unknown spot backend {name!r}; known backends: {', '.join(list_backends())} "
            f"(spots.backend selects one; a package adds one through the {ENTRY_POINT_GROUP!r} "
            f"entry-point group)") from None


def list_backends() -> list[str]:
    """Sorted names of every available backend (entry points discovered on first use)."""
    discover_backends()
    return sorted(_REGISTRY)


def discover_backends(force: bool = False) -> None:
    """Load the ``germquant.spot_backends`` entry points once (``force=True`` scans again). Every
    failure (a broken package metadata, an import error, a non-callable) is logged and skipped; an
    entry point whose name is already registered is ignored, so a built-in is never replaced."""
    global _discovered
    if _discovered and not force:
        return
    _discovered = True
    try:
        eps = list(_iter_entry_points())
    except Exception as e:  # noqa: BLE001 - discovery is an enrichment, the built-ins must still work
        log.warning("spot backend entry-point discovery failed (%s: %s); built-ins only.", type(e).__name__, e)
        return
    for ep in eps:
        name = getattr(ep, "name", None)
        value = getattr(ep, "value", "?")
        if not isinstance(name, str) or not _NAME_RE.fullmatch(name):
            log.warning("spot backend entry point with an unusable name %r (%s) ignored.", name, value)
            continue
        if name in _REGISTRY:
            log.warning("spot backend entry point %r (%s) would replace a registered backend; ignored.",
                        name, value)
            continue
        try:
            fn = ep.load()
        except Exception as e:  # noqa: BLE001 - a broken plugin must not break the pipeline
            log.warning("spot backend entry point %r (%s) failed to load (%s: %s); ignored.",
                        name, value, type(e).__name__, e)
            continue
        if not callable(fn):
            log.warning("spot backend entry point %r (%s) is not callable; ignored.", name, value)
            continue
        _REGISTRY[name] = fn
        log.info("spot backend %r loaded from entry point %s", name, value)


def _iter_entry_points() -> Iterable[Any]:
    """The raw entry points of the group (a seam the tests replace)."""
    from importlib.metadata import entry_points

    return entry_points(group=ENTRY_POINT_GROUP)


def conform_result(name: str, result: Any) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Check a backend's return value against the table contract. Missing declared columns raise;
    the declared columns are kept first (extras, if any, after them, as the append-only schema rule
    allows); a frame that already matches the contract is returned untouched."""
    try:
        per_spot, per_nuc = result
    except (TypeError, ValueError):
        raise TypeError(f"spot backend {name!r} must return (per_spot, per_nucleus) DataFrames, "
                        f"got {type(result).__name__}") from None
    per_spot = _conform(name, per_spot, PER_SPOT_COLS, "per-spot")
    per_nuc = _conform(name, per_nuc, PER_NUC_COLS, "per-nucleus")
    for what, df in (("per-spot", per_spot), ("per-nucleus", per_nuc)):
        if len(df) and (df["detector"].astype(str) != name).any():
            log.warning("spot backend %r wrote %s detector values %s; the contract is the backend name.",
                        name, what, sorted(df["detector"].astype(str).unique().tolist()))
    return per_spot, per_nuc


def _conform(name: str, df: Any, cols: list[str], what: str) -> pd.DataFrame:
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"spot backend {name!r} {what} table must be a DataFrame, got {type(df).__name__}")
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"spot backend {name!r} {what} table lacks the columns {missing}; "
                         f"the contract is {cols}")
    if list(df.columns) == cols:
        return df
    extras = [c for c in df.columns if c not in cols]
    return df[cols + extras]


def _check_name(name: Any) -> None:
    if not isinstance(name, str) or not _NAME_RE.fullmatch(name):
        raise ValueError(f"spot backend name must be a non-empty string of letters, digits, '_', '.' or '-', "
                         f"got {name!r}")


# ---- built-ins (lazy imports: detect.py imports this module, so the bodies live in sibling modules) ----

def _spotmax_backend(image, labels, spacing, *, cfg_section: dict, marker: str):
    """The validated default: the SpotMAX chain exactly as detect_spots always ran it. The resolved
    SpotMAX parameters are read from the section (detect_spots overlays its explicit arguments, so the
    code defaults are the ones detect_spots declares)."""
    from .detect import SPOTMAX_KEYS, detect_spots_spotmax

    kw = {k: cfg_section[k] for k in SPOTMAX_KEYS if k in cfg_section}
    return detect_spots_spotmax(image, labels, spacing, marker=marker, **kw)


def _log_backend(image, labels, spacing, *, cfg_section: dict, marker: str):
    """EXPERIMENTAL classical baseline (Laplacian of Gaussian); never the default."""
    from .log_detect import detect_spots_log, log_kwargs

    return detect_spots_log(image, labels, spacing, marker=marker, **log_kwargs(cfg_section))


register_backend(DEFAULT_BACKEND, _spotmax_backend)
register_backend("log", _log_backend)
