"""The public surface of germquant for GUIs, notebooks and other packages.

Stability promise. The names in ``__all__`` are the supported way to drive germquant from code: they keep
their signatures and meaning across releases of the same major version (within the 0.x line for now).
New keyword arguments and new names may be added; nothing listed here is renamed, removed or given a
different meaning without a deprecation period of at least one minor release. ``run`` is the function
``germquant run`` itself calls, so a run started from a GUI and one started from the terminal with the same
settings write the same files, and the numbers they write are gated by tests/test_golden.py. Everything
else in the package is internal and may change between releases.

Typical use::

    from germquant import api

    meta = api.read_image_metadata("gonad.ome.tif")    # channel names and voxel size, no pixels loaded
    res = api.run("gonad.ome.tif", profile="rad51_foci", out="results/gonad")
    print(res["n_nuclei"], res["qc_pass"], res["out_dir"])
    nuclei = res["tables"]["nuclei"]                   # the per-nucleus table as a DataFrame
"""
from __future__ import annotations

from pathlib import Path

from . import __version__
from .batch import collect, discover_files, enabled_stage_names, is_done, load_exclusions
from .config import Config, default_config_path, find_profile, load_config
from .io import (
    SUPPORTED_SUFFIXES,
    Stack,
    UnsupportedImageError,
    image_format,
    image_stem,
    read_image,
    read_image_metadata,
    read_nd2_metadata,
    read_stack,
)
from .models_registry import get_model
from .pipeline import process_image, resolve_model_path, spot_instances
from .provenance import read_done_marker
from .segment import segment_nuclei
from .spots import detect_spots
from .stages import STAGES, apply_cli_switches, cli_switches, stage_enabled

# one image (run, run_loaded, resolve_config, default_out_dir, process_image); configuration (Config,
# load_config, find_profile, default_config_path); reading (Stack, SUPPORTED_SUFFIXES, UnsupportedImageError,
# image_format, image_stem, read_image, read_image_metadata, read_nd2_metadata, read_stack); stages
# (STAGES, stage_enabled, cli_switches, apply_cli_switches, spot_instances); the heavy building blocks
# (segment_nuclei, detect_spots); models (resolve_model_path, get_model); batch helpers and the completion
# marker (discover_files, load_exclusions, collect, is_done, enabled_stage_names, read_done_marker).
__all__ = [
    "STAGES",
    "SUPPORTED_SUFFIXES",
    "Config",
    "Stack",
    "UnsupportedImageError",
    "__version__",
    "apply_cli_switches",
    "cli_switches",
    "collect",
    "default_config_path",
    "default_out_dir",
    "detect_spots",
    "discover_files",
    "enabled_stage_names",
    "find_profile",
    "get_model",
    "image_format",
    "image_stem",
    "is_done",
    "load_config",
    "load_exclusions",
    "process_image",
    "read_done_marker",
    "read_image",
    "read_image_metadata",
    "read_nd2_metadata",
    "read_stack",
    "resolve_config",
    "resolve_model_path",
    "run",
    "run_loaded",
    "segment_nuclei",
    "spot_instances",
    "stage_enabled",
]


def resolve_config(config: str | Path | Config | None = None, profile: str | None = None) -> Config:
    """The Config a run uses. A loaded ``Config`` passes through; a path is loaded with load_config (it may
    use ``extends:``); a ``profile`` name is looked up as ``config/profiles/<name>.yaml`` (the checkout, the
    working directory, then ``$GERMQUANT_CONFIG_ROOT``); with neither, the checkout's ``config/config.yaml``
    is used, as the drag-and-drop scripts do. Giving both raises ValueError."""
    if config is not None and profile is not None:
        raise ValueError("give either config or profile, not both")
    if isinstance(config, Config):
        return config
    if profile is not None:
        return load_config(find_profile(profile))
    if config is not None:
        return load_config(config)
    return load_config(default_config_path())


def default_out_dir(image: str | Path) -> Path:
    """Where ``run`` writes when no ``out`` is given: ``<image stem>_results`` next to the image, the
    convention of quantify.bat and quantify.command."""
    p = Path(image)
    return p.parent / f"{image_stem(p)}_results"


def run(
    image: str | Path,
    config: str | Path | Config | None = None,
    profile: str | None = None,
    out: str | Path | None = None,
    xy_stride: int | None = None,
    z_range: tuple[int, int] | list[int] | None = None,
    overrides: dict | None = None,
) -> dict:
    """Process one image exactly as ``germquant run`` does and return its result record.

    ``image`` is the path to an ``.nd2`` or, with the ``readers`` extra installed, a ``.tif`` /
    ``.ome.tif``, ``.czi``, ``.lif`` or ``.zarr`` / ``.ome.zarr`` file. ``config`` is a config YAML path or
    an already loaded Config, ``profile`` a profile name instead of it, and with neither the checkout's
    ``config/config.yaml`` is used (see resolve_config). ``out`` is the results folder (created), by
    default ``<image stem>_results`` next to the image. ``xy_stride`` downsamples x and y for a quick test
    (None or 1 = full resolution) and ``z_range`` keeps planes ``(z0, z1)``. ``overrides`` maps dotted config
    keys to values, for example ``{"spots.enabled": False}``; each is applied with ``Config.set``, so it
    enters the config hash and the provenance manifest exactly like a ``--no-<stage>`` switch on the
    command line (a Config passed in is modified in place).

    Returns the dict process_image returns: ``image_id``, ``n_nuclei``, ``qc_pass``, ``qc_flags``,
    ``out_dir``, ``tables`` (one DataFrame per output table) and ``stages`` (the per-stage record). A QC
    flag is not an error: only a failed read, segmentation or write raises.
    """
    cfg = resolve_config(config, profile)
    for key, value in (overrides or {}).items():
        cfg.set(str(key), value)
    return run_loaded(image, cfg, out, xy_stride=xy_stride, z_range=z_range)


def run_loaded(
    image: str | Path,
    cfg: Config,
    out: str | Path | None = None,
    *,
    xy_stride: int | None = None,
    z_range: tuple[int, int] | list[int] | None = None,
) -> dict:
    """``run`` once the config is resolved and every override applied: the single call that the command
    line and the API share, so there is one code path to a results folder."""
    out_dir = Path(out) if out is not None else default_out_dir(image)
    stride = int(xy_stride) if xy_stride else 1
    zr = (int(z_range[0]), int(z_range[1])) if z_range else None
    # process_image writes run_manifest.json itself (run geometry and per-stage hashes), hence prov=None
    return process_image(image, cfg, out_dir, xy_stride=stride, z_range=zr, prov=None)
