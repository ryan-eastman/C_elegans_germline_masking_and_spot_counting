"""Pluggable nucleus segmentation backends (docs/BACKENDS.md).

A backend is any callable with this signature:

    backend(dna, spacing, *, model, diameter_um, device, cfg_section) -> labels

* ``dna``: the DNA (DAPI) channel as a (Z, Y, X) array, as read from the .nd2 (uint16 or float32).
* ``spacing``: (dz, dy, dx) voxel size in micrometres from the .nd2 metadata.
* ``model``: the resolved value of ``segmentation.nuclei.cellpose_model`` (a model file path, resolved
  against the repository or ``$GERMQUANT_MODELS``, or a bare built-in name such as 'cpsam'). A backend
  may use it as its own model reference or ignore it.
* ``diameter_um``: ``segmentation.nuclei.diameter_um``, the expected nucleus diameter.
* ``device``: ``segmentation.nuclei.device`` ('cuda' | 'mps' | 'cpu' | 'auto' | None). None means
  auto; ``$GERMQUANT_DEVICE`` overrides it (``germquant.device.select_device``).
* ``cfg_section``: a copy of the whole ``segmentation.nuclei`` config section, so a backend can read
  keys of its own without a change to this package.
* returns an integer label image of shape (Z, Y, X): 0 is background, every nucleus one label.

``segment_nuclei`` applies the minimum-volume filter and the consecutive relabel afterwards, exactly as
it always has for the Cellpose path, so a backend only produces instance labels. Built-ins: 'cellpose'
(the Cellpose-SAM path with the explicit device and the CPU retry) and 'classical' (the spacing-aware
watershed). An installed package adds a backend by declaring an entry point in the group
``germquant.segment_backends``; the entry point's name is the backend name. A backend that fails to
load is logged and skipped, never fatal, so a broken plugin cannot stop a run that does not use it.
"""
from __future__ import annotations

import logging
from collections.abc import Callable
from importlib import metadata

import numpy as np

log = logging.getLogger(__name__)

ENTRY_POINT_GROUP = "germquant.segment_backends"
DEFAULT_BACKEND = "cellpose"

Backend = Callable[..., np.ndarray]


class UnknownBackendError(LookupError):
    """A backend name that is neither built in, registered in this process, nor provided by an
    installed package's entry point."""


# ---- the two built-ins ----------------------------------------------------------------------

def cellpose_backend(dna, spacing, *, model, diameter_um, device, cfg_section) -> np.ndarray:
    """The Cellpose-SAM path: ``nuclei._cellpose`` with the explicit device from ``germquant.device``
    (honours ``$GERMQUANT_DEVICE``) and one retry on the CPU when the accelerator fails (out of memory,
    an operator Metal cannot run). ``cfg_section`` is unused: every Cellpose knob is an argument."""
    from ..device import select_device
    from . import nuclei  # looked up at call time, so a test can monkeypatch nuclei._cellpose

    _ = cfg_section
    try:
        return nuclei._cellpose(dna, spacing, model, diameter_um, device_pref=device)
    except Exception as e:  # accelerator trouble: the same model on the CPU first
        if select_device(device) == "cpu":
            raise
        log.warning("Cellpose failed on %s (%s: %s); retrying on the CPU.", select_device(device),
                    type(e).__name__, e)
        return nuclei._cellpose(dna, spacing, model, diameter_um, device_pref="cpu")


def classical_backend(dna, spacing, *, model, diameter_um, device, cfg_section) -> np.ndarray:
    """The spacing-aware watershed (``nuclei._classical``): a baseline that needs no GPU and no model.
    ``model``, ``device`` and ``cfg_section`` are unused."""
    from . import nuclei

    _ = (model, device, cfg_section)
    return nuclei._classical(dna, spacing, diameter_um)


_BUILTINS: dict[str, Backend] = {"cellpose": cellpose_backend, "classical": classical_backend}
_REGISTRY: dict[str, Backend] = dict(_BUILTINS)
_SOURCES: dict[str, str] = {name: "built-in" for name in _BUILTINS}
_FAILED: dict[str, str] = {}      # entry points that could not be loaded, name -> reason (for error messages)
_STATE: dict[str, bool] = {"discovered": False}      # entry points are scanned once per process


# ---- the registry ---------------------------------------------------------------------------

def register_backend(name: str, fn: Backend, *, replace: bool = False) -> None:
    """Make ``fn`` available as backend ``name`` in this process (an installed package does this
    through its entry point; call it directly from a script, a notebook or a test). Registering a name
    twice is an error unless ``replace=True``, which also lets a test stand in for a built-in."""
    if not isinstance(name, str):
        raise TypeError(f"a backend name must be a string, got {type(name).__name__}")
    name = name.strip()
    if not name:
        raise ValueError("a backend name must not be empty")
    if not callable(fn):
        raise TypeError(f"backend {name!r} must be a callable, got {type(fn).__name__}")
    if name in _REGISTRY and not replace:
        raise ValueError(f"backend {name!r} is already registered ({_SOURCES.get(name, 'registered')}); "
                         f"pass replace=True to override it")
    _REGISTRY[name] = fn
    _SOURCES[name] = "registered"


def unregister_backend(name: str) -> None:
    """Undo ``register_backend``. A built-in name goes back to its built-in implementation; it is
    never removed."""
    if name in _BUILTINS:
        _REGISTRY[name] = _BUILTINS[name]
        _SOURCES[name] = "built-in"
        return
    if name not in _REGISTRY:
        raise UnknownBackendError(_unknown_message(name))
    del _REGISTRY[name]
    del _SOURCES[name]


def get_backend(name: str | None = None) -> Backend:
    """The callable registered as ``name`` (None or '' means the default, 'cellpose'): a built-in, a
    backend registered in this process, or an installed package's entry point (discovered on the first
    call). An unknown name raises ``UnknownBackendError`` naming what is available; it is never
    silently replaced by another backend."""
    discover_entry_points()
    key = _normalise(name)
    try:
        return _REGISTRY[key]
    except KeyError:
        raise UnknownBackendError(_unknown_message(key)) from None


def list_backends() -> list[str]:
    """Sorted names of every available backend (built-ins, registered, installed entry points)."""
    discover_entry_points()
    return sorted(_REGISTRY)


def backend_sources() -> dict[str, str]:
    """``{name: origin}`` for every available backend: 'built-in', 'registered' (register_backend in
    this process) or 'entry point module:attr' (an installed package)."""
    discover_entry_points()
    return {name: _SOURCES[name] for name in sorted(_REGISTRY)}


def discover_entry_points(*, force: bool = False) -> dict[str, str]:
    """Register the backends that installed packages declare under ``ENTRY_POINT_GROUP``.

    Runs once per process, lazily on the first lookup; ``force=True`` scans again (for a package
    installed after the first lookup). Returns ``{backend name: 'module:attr'}`` for the entry points
    this call loaded. A load failure (import error, missing attribute, not a callable) is logged with
    the entry point's origin and skipped, never raised, so a broken plugin cannot stop a run that does
    not use it. A name that is already registered (a built-in, or one registered earlier in this
    process) keeps its implementation and the entry point is skipped with a warning.
    """
    if _STATE["discovered"] and not force:
        return {}
    _STATE["discovered"] = True
    loaded: dict[str, str] = {}
    try:
        eps = _iter_entry_points()
    except Exception as e:  # noqa: BLE001 - a damaged installation must not stop the run
        log.warning("could not read the %r entry points (%s: %s); only the built-in segmentation "
                    "backends are available", ENTRY_POINT_GROUP, type(e).__name__, e)
        return loaded
    for ep in eps:
        name = str(ep.name).strip()
        origin = f"entry point {ep.value}"
        if name in _REGISTRY:
            if _SOURCES.get(name) != origin:          # a forced rescan meets its own earlier loads
                log.warning("segmentation backend %r (%s) is already registered as %s; keeping that one",
                            name, ep.value, _SOURCES.get(name))
            continue
        try:
            fn = ep.load()
        except Exception as e:  # noqa: BLE001 - never fatal: logged, and that backend is unavailable
            log.warning("segmentation backend %r (%s) failed to load and is unavailable: %s: %s",
                        name, ep.value, type(e).__name__, e)
            _FAILED[name] = f"{ep.value}: {type(e).__name__}: {e}"
            continue
        if not callable(fn):
            log.warning("segmentation backend %r (%s) is not callable and is unavailable", name, ep.value)
            _FAILED[name] = f"{ep.value}: not callable"
            continue
        _REGISTRY[name] = fn
        _SOURCES[name] = origin
        loaded[name] = ep.value
    return loaded


def _iter_entry_points() -> list:
    """The installed entry points of our group (a seam for tests)."""
    return list(metadata.entry_points(group=ENTRY_POINT_GROUP))


def _normalise(name) -> str:
    if name is None:
        return DEFAULT_BACKEND
    key = str(name).strip()
    return key or DEFAULT_BACKEND


def _unknown_message(name) -> str:
    failed = ""
    if _FAILED:
        failed = (" Installed backends that failed to load: "
                  + "; ".join(f"{k} ({v})" for k, v in sorted(_FAILED.items())) + ".")
    return (f"unknown nucleus segmentation backend {name!r}; available: {', '.join(sorted(_REGISTRY))}.{failed} "
            f"Check the spelling of segmentation.nuclei.backend; a third-party backend must be installed "
            f"in this Python environment and declare an entry point in the group {ENTRY_POINT_GROUP!r} "
            f"(docs/BACKENDS.md).")
