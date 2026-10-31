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
from .batch import collect, discover_files, discover_options, enabled_stage_names, is_done, load_exclusions
from .config import Config, default_config_path, find_profile, load_config
from .import_labels import import_labels
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
from .io.sample_sheet import (
    SampleRow,
    SampleSheet,
    SampleSheetError,
    SampleSheetWarning,
    read_sample_sheet,
    write_sample_sheet,
)
from .models_registry import get_model, known_models, register_local_model
from .packs import fetch_example, fetch_pack_models, find_pack, install_pack, list_packs, load_pack, pack_status
from .pipeline import process_image, resolve_model_path, spot_instances
from .project import (
    Project,
    ProjectError,
    ProjectLockedError,
    RunError,
    add_images,
    break_lock,
    collect_project,
    export_commands,
    list_rounds,
    new_project,
    new_round,
    open_project,
    preview_roles,
    project_config,
    project_status,
    revert_recipe,
    run_project,
    save_project,
    use_model,
    write_recipe,
)
from .project.headers import read_header
from .provenance import read_done_marker
from .segment import segment_nuclei
from .spots import detect_spots
from .stages import STAGES, apply_cli_switches, cli_switches, stage_enabled
from .timelapse import parse_timepoints, process_series, register_series
from .tune import GridError, TuneRefused, TuningWarning, apply_tuning, tune_spots

# one image (run, run_loaded, resolve_config, default_out_dir, process_image); a time-lapse (process_series,
# parse_timepoints: the loop `germquant run --t` runs; register_series: `germquant register`, the drift of a
# series as a shift table, docs/TIMELAPSE.md); configuration (Config,
# load_config, find_profile, default_config_path); reading (Stack, SUPPORTED_SUFFIXES, UnsupportedImageError,
# image_format, image_stem, read_image, read_image_metadata, read_nd2_metadata, read_stack); stages
# (STAGES, stage_enabled, cli_switches, apply_cli_switches, spot_instances); the heavy building blocks
# (segment_nuclei, detect_spots); models (resolve_model_path, get_model); batch helpers and the completion
# marker (discover_files, load_exclusions, collect, is_done, enabled_stage_names, read_done_marker); projects
# (germquant.project, docs/PROJECTS.md: Project, new_project to run_project and export_commands, the light header
# read_header; importing germquant.project alone loads no torch); spot tuning from truth points (tune_spots and
# apply_tuning: `germquant tune-spots`, docs/TUNING.md, with TuneRefused, GridError and TuningWarning; germquant.tune
# holds the tuners, grids and point files); training pairs from labels drawn in another program (import_labels:
# `germquant import-labels`, docs/TRAINING.md); the sample sheet (SampleRow, SampleSheet, read_sample_sheet,
# write_sample_sheet, with SampleSheetError and SampleSheetWarning: `--samples`, docs/PROJECTS.md "Sample sheet";
# run and run_loaded take a row as sample=); batch formats (discover_options: the io.input_suffixes keywords
# `germquant batch` gives discover_files).
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
    "parse_timepoints",
    "process_image",
    "process_series",
    "read_done_marker",
    "read_image",
    "read_image_metadata",
    "read_nd2_metadata",
    "read_stack",
    "register_series",
    "resolve_config",
    "resolve_model_path",
    "run",
    "run_loaded",
    "segment_nuclei",
    "spot_instances",
    "stage_enabled",
    # system packs (docs/PACKS.md); germquant.packs offers the same listing and checks without loading torch
    "list_packs",
    "find_pack",
    "load_pack",
    "pack_status",
    "fetch_pack_models",
    "fetch_example",
    "install_pack",
    "list_segment_backends",
    "list_spot_backends",
    "known_models",
    "register_local_model",
    # projects (docs/PROJECTS.md); germquant.project offers the same names, and its light entry the writes, without torch
    "Project",
    "ProjectError",
    "ProjectLockedError",
    "RunError",
    "new_project",
    "open_project",
    "save_project",
    "add_images",
    "read_header",
    "project_config",
    "preview_roles",
    "write_recipe",
    "revert_recipe",
    "use_model",
    "run_project",
    "project_status",
    "collect_project",
    "export_commands",
    "break_lock",
    "list_rounds",
    "new_round",
    # spot tuning (docs/TUNING.md); importing germquant.tune alone loads no torch, and each tuner imports the pipeline;
    # TuneRefused (a ValueError) and GridError (a TuneRefused) are the refusals of tune_spots and apply_tuning, and
    # TuningWarning the warning of apply_tuning (PLAN.md 9.4 #12)
    "tune_spots",
    "apply_tuning",
    "TuneRefused",
    "GridError",
    "TuningWarning",
    # training pairs (docs/TRAINING.md); importing germquant.import_labels alone loads no torch, and the call
    # imports the pipeline before it reads the image
    "import_labels",
    # the sample sheet (docs/PROJECTS.md, "Sample sheet"); importing germquant.io.sample_sheet alone loads no torch
    "SampleRow",
    "SampleSheet",
    "read_sample_sheet",
    "write_sample_sheet",
    # SampleSheetError (a ValueError): a sheet read_sample_sheet cannot read or write_sample_sheet refuses, or a row
    # given to run with t; SampleSheetWarning (a UserWarning): each warning of a sheet that is read (PLAN.md 9.4 #13)
    "SampleSheetError",
    "SampleSheetWarning",
    # batch formats (README, "Other file formats"; PLAN.md 9.1 with 9.4 #13): {} for a recipe without
    # io.input_suffixes, so discover_files keeps the .nd2 files it always has
    "discover_options",
]


def list_segment_backends() -> list[str]:
    """The segmentation backend names this environment declares (built-ins, in-process registrations and
    installed entry points), read from metadata without loading any backend."""
    from .segment.backends import declared_backend_names

    return declared_backend_names()


def list_spot_backends() -> list[str]:
    """The spot backend names this environment declares, read from metadata without loading any backend."""
    from .spots.backends import declared_backend_names

    return declared_backend_names()


def resolve_config(config: str | Path | Config | None = None, profile: str | None = None,
                   pack: str | None = None, recipe: str | None = None) -> Config:
    """The Config a run uses, from exactly one source. A loaded ``Config`` passes through; a path is loaded
    with load_config (it may use ``extends:``); a ``profile`` name is looked up as
    ``config/profiles/<name>.yaml`` (the checkout, the working directory, then ``$GERMQUANT_CONFIG_ROOT``); a
    ``pack`` (``id`` or ``id@version``) loads one of its recipes (``recipe``, else the pack's default) through
    load_config and sets ``cfg.pack_record`` (docs/PACKS.md): the config_hash equals that of ``config=<the
    recipe file>``. With no source, the checkout's ``config/config.yaml`` is used, as the drag-and-drop scripts
    do. More than one source, or ``recipe`` without ``pack``, raises ValueError."""
    given = [name for name, v in (("config", config), ("profile", profile), ("pack", pack)) if v is not None]
    if len(given) > 1:
        raise ValueError(f"give either config, profile or pack, not both ({' and '.join(given)} given)")
    if recipe is not None and pack is None:
        raise ValueError("recipe needs pack")
    if pack is not None:
        from .packs import find_pack as _find_pack
        from .packs import pack_record, recipe_profile_path

        p = _find_pack(pack)
        r = p.recipe(recipe)
        cfg = load_config(recipe_profile_path(p, r))
        cfg.pack_record = pack_record(p, r.id)
        return cfg
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
    t: int | None = None,
    pack: str | None = None,
    recipe: str | None = None,
    sample=None,
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
    command line (a Config passed in is modified in place). ``t`` processes one timepoint of a time series:
    the stages run on that timepoint's stack and the image id becomes ``<stem>_t<NNN>`` (docs/TIMELAPSE.md);
    None, the default, for a fixed image. ``pack`` (``id`` or ``id@version``) and ``recipe`` run a system
    pack's recipe instead of config or profile (docs/PACKS.md), as ``germquant run --pack`` does; the run
    also records the pack and the runtime in its manifest. This call does not run the pack check: call
    ``pack_status`` first, as the command line does. ``sample`` is the image's row of a sample sheet
    (``germquant.io.sample_sheet.read_sample_sheet(path).row(image_id, image)``, a SampleRow; docs/PROJECTS.md,
    "Sample sheet"), as ``germquant run --samples`` passes it: its design columns replace the file name's; None, the
    default, for none. A sample row with ``t`` is refused.

    Returns the dict process_image returns: ``image_id``, ``n_nuclei``, ``qc_pass``, ``qc_flags``,
    ``out_dir``, ``tables`` (one DataFrame per output table) and ``stages`` (the per-stage record). A QC
    flag is not an error: only a failed read, segmentation or write raises.
    """
    cfg = resolve_config(config, profile, pack=pack, recipe=recipe)
    for key, value in (overrides or {}).items():
        cfg.set(str(key), value)
    return run_loaded(image, cfg, out, xy_stride=xy_stride, z_range=z_range, t=t,
                      **({} if sample is None else {"sample": sample}))


def run_loaded(
    image: str | Path,
    cfg: Config,
    out: str | Path | None = None,
    *,
    xy_stride: int | None = None,
    z_range: tuple[int, int] | list[int] | None = None,
    t: int | None = None,
    sample=None,
) -> dict:
    """``run`` once the config is resolved and every override applied: the single call that the command
    line and the API share, so there is one code path to a results folder. ``sample``: the image's sample sheet
    row (see ``run``)."""
    out_dir = Path(out) if out is not None else default_out_dir(image)
    stride = int(xy_stride) if xy_stride else 1
    zr = (int(z_range[0]), int(z_range[1])) if z_range else None
    # process_image writes run_manifest.json itself (run geometry and per-stage hashes), hence prov=None;
    # t goes along only when given, so a fixed image is the call it has always been
    kw = {} if t is None else {"t": t}
    if sample is not None:
        kw["sample"] = sample              # a sample sheet row, passed on only when given (D40)
    return process_image(image, cfg, out_dir, xy_stride=stride, z_range=zr, prov=None, **kw)
