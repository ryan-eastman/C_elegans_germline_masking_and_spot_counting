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
the pre-registry code), ``log`` (an EXPERIMENTAL classical Laplacian-of-Gaussian baseline that is
never selected unless ``spots.backend: log`` is set) and ``net3d`` (an EXPERIMENTAL learned 3D U-Net
heatmap detector, ``spots.backend: net3d`` with ``spots.net3d.model``; it needs torch).

Two optional attributes of a backend callable extend the contract without touching the others (an
append-only change: a backend without them is called exactly as before):

* ``input_roles``: channel roles beyond the foci channel the backend wants (net3d: ``("dna",)``). The
  caller passes those the image has as ``detect_spots(..., channels={role: array})`` and the backend
  receives them as the keyword ``channels``; a backend without the attribute never sees that keyword.
* ``provenance(cfg_section) -> dict``: what a run records before any pixel is read (net3d: its model
  file's sha256), which the pipeline adds to the spots stage sub-hash and the run manifest
  (``spot_backend``) only when the dict is not empty, so the default run's hashes do not move.

Third-party packages add backends through the
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


def declared_backend_names() -> list[str]:
    """Sorted names of every spot backend this environment declares, WITHOUT loading any: the built-ins,
    the in-process registrations and the entry-point names read from the installed packages' metadata
    (no module is imported, so nothing pulls in SpotMAX, CuPy or torch). pack status checks a pack's
    needed backends with it (docs/DESIGN_system_packs.md, section 5.1)."""
    names = set(_REGISTRY)
    try:
        names |= {str(getattr(ep, "name", "")).strip() for ep in _iter_entry_points()}
    except Exception as e:  # noqa: BLE001 - a damaged installation must not stop a listing
        log.warning("could not read the %r entry points (%s: %s)", ENTRY_POINT_GROUP, type(e).__name__, e)
    return sorted(n for n in names if n and _NAME_RE.fullmatch(n))


def backend_input_roles(name: str) -> tuple[str, ...]:
    """The channel roles beyond the foci channel that backend ``name`` asks for, from its optional
    ``input_roles`` attribute (net3d: ``("dna",)``); empty for spotmax and log. The caller passes the
    roles the image has as ``detect_spots(..., channels={role: array})``. ValueError for an unknown
    name, as ``get_backend``."""
    roles = getattr(get_backend(name), "input_roles", ()) or ()
    return tuple(str(r) for r in roles)


def backend_provenance(name: str, cfg_section: dict | None) -> dict:
    """What a run records about backend ``name`` before any pixel is read, from its optional
    ``provenance(cfg_section)`` attribute (net3d: the model file and its sha256). {} for a backend
    without one and for an unknown name (the spots stage reports that one when it runs); a failing hook
    is logged and recorded as ``provenance_error``, never fatal."""
    try:
        fn = get_backend(name)
    except (TypeError, ValueError):
        return {}
    hook = getattr(fn, "provenance", None)
    if hook is None:
        return {}
    try:
        rec = hook(dict(cfg_section or {}))
    except Exception as e:  # noqa: BLE001 - provenance is a record, it must not stop a run
        log.warning("spot backend %r: provenance hook failed (%s: %s)", name, type(e).__name__, e)
        return {"backend": name, "provenance_error": f"{type(e).__name__}: {e}"}
    return dict(rec or {})


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


def _net3d_backend(image, labels, spacing, *, cfg_section: dict, marker: str, channels: dict | None = None):
    """EXPERIMENTAL learned detector (a trained net3d 3D U-Net, spots.net3d.model); never the default.
    Needs torch; reads the DNA channel from ``channels`` when the model was trained with it."""
    from .net3d.backend import detect_spots_net3d, net3d_kwargs

    return detect_spots_net3d(image, labels, spacing, marker=marker, channels=channels, **net3d_kwargs(cfg_section))


def _net3d_provenance(cfg_section: dict) -> dict:
    from .net3d.backend import net3d_provenance

    return net3d_provenance(cfg_section)


_net3d_backend.input_roles = ("dna",)           # passed when the image has it; a foci-only model ignores it
_net3d_backend.provenance = _net3d_provenance   # the model file's sha256 enters the spots sub-hash

register_backend(DEFAULT_BACKEND, _spotmax_backend)
register_backend("log", _log_backend)
register_backend("net3d", _net3d_backend)
