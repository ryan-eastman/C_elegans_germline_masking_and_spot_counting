"""Whether a pack recipe can run here: `germquant pack check` (docs/DESIGN_system_packs.md, section 3.9).

The shallow check imports nothing heavy (no torch, Cellpose, SpotMAX or CuPy):

* the engine version lies inside requires.germquant;
* the pack's files match its SHA256SUMS, and a root: checkout pack's engine files match their digests;
* the recipe loads, its config hash equals expected_config_hash, and each of its tunables is a key it sets,
  with a scalar value (packs.validate_pack);
* each needed backend name is declared (the built-ins plus entry-point names read from metadata), and each
  needed module is found by importlib.util.find_spec (the built-in cellpose name is always declared, so the
  module check is what greys a pack in an environment without Cellpose or torch);
* each model is present with the right size; its sha256 is verified when its size or mtime differs from its
  verification record, and trusted otherwise;
* the readers extra is present when the input is not an .nd2;
* with an image given, its pixel size and z step fall inside ``expects`` (a warning, never a block). The
  shallow check reads only an .nd2 header (the nd2 package) or a TIFF header (tifffile,
  io.ome_header): bioio would import dask.array, and dask imports CuPy whenever CuPy is installed, before
  torch. A CZI, LIF or OME-Zarr header is left to the deep check, with a warning saying so.

``deep`` imports germquant.pipeline first (torch before CuPy), then loads the backends, reports the device
that will run, with a warning when the recipe recommends an accelerator and none is present, and reads the
image header the shallow check left.
"""
from __future__ import annotations

import importlib.util
import logging
from pathlib import Path

from .format import Pack, lf_sha256, verify_sums

log = logging.getLogger(__name__)

# the extra (pyproject optional-dependencies) that brings each module a recipe may need
MODULE_EXTRA = {"torch": "gpu", "cellpose": "gpu", "spotmax": "spots", "cellacdc": "spots", "skan": "sc",
                "laptrack": "timelapse", "bioio": "readers", "h5py": "imaris", "snakemake": "workflow"}


def _module_present(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def _install_fix(module: str) -> str:
    extra = MODULE_EXTRA.get(module)
    if extra:
        return (f"install the {extra!r} extra into this environment (from the engine checkout: "
                f"uv pip install -e \".[{extra}]\")")
    return f"install the Python package that provides {module!r} into this environment"


def engine_files_problems(pack: Pack) -> list[str]:
    """A root: checkout pack records the LF-normalised sha256 of each engine file it references; problems
    when one is missing or differs (the frozen engine files never change, rule R1)."""
    if pack.root != "checkout":
        return []
    from .discovery import engine_root_for

    try:
        root = engine_root_for(pack)
    except FileNotFoundError as e:
        return [str(e)]
    out = []
    for rel, want in sorted((pack.data.get("engine_files") or {}).items()):
        f = root / rel
        if not f.is_file():
            out.append(f"engine file {rel} is missing under {root}")
        elif lf_sha256(f) != str(want):
            out.append(f"engine file {rel} differs from the digest the pack records")
    return out


def _read_header(image, *, deep: bool = False) -> dict | None:
    """The image header _expects_warnings compares, or None when only the deep check may read it. An .nd2 goes
    through io.read_image_metadata (read_nd2_metadata: the nd2 package reads the header without dask or bioio)
    and a TIFF through io.ome_header.read_tiff_header (tifffile only). bioio imports dask.array, and dask
    imports CuPy whenever CuPy is installed, before torch (rule R5, gate G6), so another format is read only
    with ``deep``, which the deep check passes after it has imported the pipeline (torch first)."""
    from ..io.formats import image_format
    from ..io.ome_header import is_tiff, read_tiff_header

    if deep or image_format(image) == "nd2":
        from ..io import read_image_metadata

        return read_image_metadata(image)
    if is_tiff(image):
        return read_tiff_header(image)
    return None


def _expects_warnings(pack: Pack, image, *, deep: bool = False) -> tuple[list[dict], dict | None]:
    """``(warnings, pending)``: the image's pixel size and z step against the pack's ``expects`` (warnings
    only; nothing is read when the pack expects nothing). ``pending`` is the warning that says the header was
    left to the deep check (see _read_header), else None."""
    from ..io.formats import image_suffix

    exp = pack.data.get("expects") or {}
    exp = exp if isinstance(exp, dict) else {}
    checks = [(key, what) for key, what in (("pixel_um", "pixel size"), ("z_step_um", "z step")) if exp.get(key)]
    if not checks:
        return [], None
    try:
        meta = _read_header(image, deep=deep)
    except Exception as e:  # noqa: BLE001 - the header is advisory here
        return [{"what": f"could not read the header of {image}: {type(e).__name__}: {e}",
                 "fix": "check the file and the readers extra"}], None
    if meta is None:
        pending = {"what": f"the pixel size and z step of {Path(image).name} were not checked: without --deep, pack "
                           f"check reads .nd2 and TIFF headers only (a {image_suffix(image)} header needs the readers, "
                           f"which load dask and, where installed, CuPy)",
                   "fix": f"run `germquant pack check {pack.ref} --image FILE --deep`, which reads it after loading "
                          f"the pipeline"}
        return [pending], pending
    if not meta.get("spacing_ok", True):
        return [], None
    dz, dy, dx = (float(v) for v in meta.get("spacing", (float("nan"),) * 3))
    values = {"pixel_um": dy, "z_step_um": dz}
    out = []
    for key, what in checks:
        lo, hi = float(exp[key][0]), float(exp[key][1])
        val = values[key]
        if not lo <= val <= hi:
            out.append({"what": f"the image's {what} {val:.4g} um lies outside the pack's expected range "
                                f"[{lo:g}, {hi:g}] um", "fix": "check that this pack suits the image (a warning only)"})
    return out, None


def pack_status(pack, recipe: str | None = None, image=None, deep: bool = False) -> dict:
    """``{"pack", "version", "recipe", "usable", "blockers", "warnings", "models", "device"}``. Each blocker and
    warning is ``{"what", "fix"}``. ``pack`` is a Pack or a name (``id`` or ``id@version``)."""
    from .. import __version__
    from ..config import load_config
    from .blobs import check_blob
    from .discovery import find_pack, recipe_profile_path
    from .validate import tunable_problems

    blockers: list[dict] = []
    warnings: list[dict] = []
    p = pack if isinstance(pack, Pack) else find_pack(pack)
    out = {"pack": p.id, "version": p.version, "recipe": None, "usable": False, "blockers": blockers,
           "warnings": warnings, "models": [], "device": None, "engine": __version__}
    for prob in p.problems:
        blockers.append({"what": f"pack manifest: {prob}", "fix": "use a corrected pack version"})
    try:
        r = p.recipe(recipe)
    except FileNotFoundError as e:
        blockers.append({"what": str(e), "fix": "pick one of the pack's recipes with --recipe"})
        return out
    out["recipe"] = r.id

    spec = (p.data.get("requires") or {}).get("germquant")
    if spec:
        from .discovery import engine_compatible

        if not engine_compatible(p):
            blockers.append({"what": f"engine {__version__} is outside the pack's range {spec}",
                             "fix": "update the engine checkout to a version in that range (git pull), or pick a "
                                    "version of the pack made for this engine"})
    sums = verify_sums(p.folder)
    if sums:
        blockers.append({"what": f"the pack's files do not match its SHA256SUMS ({sums[0]}"
                                 + (f" and {len(sums) - 1} more" if len(sums) > 1 else "") + ")",
                         "fix": f"reinstall pack {p.ref} (germquant pack install SRC --replace)"})
    for prob in engine_files_problems(p):
        blockers.append({"what": prob, "fix": "use the engine checkout this pack was made for"})

    try:
        cfg = load_config(recipe_profile_path(p, r))
        if cfg.hash != r.expected_config_hash:
            blockers.append({"what": f"recipe {r.id} hashes to {cfg.hash}, the pack expects {r.expected_config_hash}",
                             "fix": "the recipe text or its channel map changed: reinstall the pack version, or ask "
                                    "its maintainer for a new version"})
    except Exception as e:  # noqa: BLE001 - any failure to load is a blocker with its message
        cfg = None
        blockers.append({"what": f"recipe {r.id} does not load: {type(e).__name__}: {e}",
                         "fix": "reinstall the pack, or check its extends chain (germquant pack verify)"})
    if cfg is not None:
        for prob in tunable_problems(r, cfg.as_dict()):             # what packs.validate_pack refuses
            blockers.append({"what": prob, "fix": "use a corrected pack version (germquant pack lint names every "
                                                  "faulty tunable)"})

    from ..segment.backends import declared_backend_names as seg_declared
    from ..spots.backends import declared_backend_names as spot_declared

    seg_have, spot_have = set(seg_declared()), set(spot_declared())
    for name in r.need("segment_backends"):
        if name not in seg_have:
            blockers.append({"what": f"segmentation backend {name!r} is not installed",
                             "fix": f"install the package that provides the {name!r} backend (an entry point in "
                                    f"germquant.segment_backends) into this environment"})
    for name in r.need("spot_backends"):
        if name not in spot_have:
            blockers.append({"what": f"spot backend {name!r} is not installed",
                             "fix": f"install the package that provides the {name!r} backend (an entry point in "
                                    f"germquant.spot_backends) into this environment"})
    for mod in r.need("modules"):
        if not _module_present(mod):
            blockers.append({"what": f"Python module {mod!r} is not installed", "fix": _install_fix(mod)})
    for mname in r.need("models"):
        mspec = p.model(mname)
        if mspec is None:
            continue                                   # a manifest problem, reported above
        st = check_blob(mspec)
        out["models"].append({"name": mname, **st})
        if st["state"] in ("missing", "wrong_size", "mismatch"):
            blockers.append({"what": f"model {mname} is not usable: {st['detail']} ({st['path']})",
                             "fix": f"germquant pack fetch-models {p.id}"
                                    + (f" --recipe {r.id}" if len(p.recipes) > 1 else "")})
    pending = None
    if image is not None:
        from ..io.formats import image_format

        try:
            fmt = image_format(image)
        except ValueError as e:
            fmt = None
            blockers.append({"what": str(e), "fix": "use a supported image format"})
        if fmt == "bioio" and not _module_present("bioio"):
            blockers.append({"what": f"{Path(image).name} is not an .nd2 and the readers extra is not installed",
                             "fix": _install_fix("bioio")})
        elif fmt is not None:
            w, pending = _expects_warnings(p, image)
            warnings.extend(w)
    if deep and not blockers:
        _deep_check(p, r, cfg, out, image=image, pending=pending)
    out["usable"] = not blockers
    return out


def _deep_check(pack: Pack, recipe, cfg, out: dict, *, image=None, pending: dict | None = None) -> None:
    """Import the pipeline first (torch before CuPy, the cuBLAS load order), load the needed backends, report
    the device that will run, and read the image header the shallow check left (``pending``)."""
    import germquant.pipeline  # noqa: F401  (side effect: torch claims cuBLAS before any backend import)

    from ..device import select_device
    from ..segment.backends import get_backend as get_seg
    from ..spots.backends import get_backend as get_spot

    if pending is not None:                        # torch is loaded now, so the readers may load CuPy after it
        out["warnings"].remove(pending)
        out["warnings"].extend(_expects_warnings(pack, image, deep=True)[0])

    for name in recipe.need("segment_backends"):
        try:
            get_seg(name)
        except Exception as e:  # noqa: BLE001 - reported, not raised
            out["blockers"].append({"what": f"segmentation backend {name!r} fails to load: {type(e).__name__}: {e}",
                                    "fix": "reinstall the package that provides it"})
    for name in recipe.need("spot_backends"):
        try:
            get_spot(name)
        except Exception as e:  # noqa: BLE001
            out["blockers"].append({"what": f"spot backend {name!r} fails to load: {type(e).__name__}: {e}",
                                    "fix": "reinstall the package that provides it"})
    prefer = cfg.get("segmentation.nuclei.device") if cfg is not None else None
    dev = select_device(prefer, quiet=True)
    out["device"] = dev
    if str(recipe.needs.get("accelerator", "optional")) == "recommended" and dev == "cpu":
        out["warnings"].append({"what": "the recipe recommends an accelerator (CUDA or Apple Metal) and none is "
                                        "available: it will run on the CPU, slowly",
                                "fix": "run `germquant check-gpu`; see docs/MAC_METAL.md or the README"})
