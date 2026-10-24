"""The model path guard: a configured Cellpose model path that names no file stops the run.

``pipeline.resolve_model_path`` hands its value to ``model_file_or_refusal`` wherever it would otherwise return
a path it did not find as a file: an absolute path, and its final fall-through. Given to Cellpose 4
(``CellposeModel(pretrained_model=...)``), a path that does not exist is taken for an unknown model name:
Cellpose logs "pretrained model ... not found, using default model" and segments with its stock weights
(cpsam_v2), while the run records model_sha256 "absent" and segmentation_method "cellpose". Here such a value
raises FileNotFoundError instead, naming the places looked in and the fixes.

Every other value comes back exactly as it was given: a file (a relative one counts from the working directory,
where Cellpose and provenance.file_sha256 look), a bare name (resolve_model_path has already refused every bare
name Cellpose would replace), and every value the built-in Cellpose backend never loads: a method the segment
stage runs as the classical watershed (classical, or any value other than auto and cellpose), another backend
(classical, nuclei3d, a plugin), an unknown backend name (the segment stage refuses it), and, under method auto,
an installation without Cellpose, where the run goes on with the classical watershed as it always did.

Only config-driven runs come through here: process_image (``germquant run``, ``api.run``), ``germquant batch``
and the time-lapse loop. A direct call of a building block (``api.segment_nuclei``, ``planar.segment_planar``,
a backend called with ``model=``) hands its model value to Cellpose unchecked, as before.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path


def cellpose_installed() -> bool:
    """Whether ``import cellpose`` would find a package, without importing it (no torch or Cellpose import, so
    the import order of a run is unchanged). A module already in sys.modules counts, a test's stand-in included;
    an import blocked with ``sys.modules["cellpose"] = None`` does not."""
    if "cellpose" in sys.modules:
        return sys.modules["cellpose"] is not None
    try:
        return importlib.util.find_spec("cellpose") is not None
    except (ImportError, ValueError):
        return False


def model_file_or_refusal(cfg, model_path):
    """``model_path`` unchanged when it is a file, a bare name, or a value the built-in Cellpose backend will not
    load (the module docstring lists them). Otherwise FileNotFoundError naming the places resolve_model_path
    looked in, in its order, what Cellpose would have done, and the fixes."""
    p = Path(str(model_path))
    if p.is_file() or (len(p.parts) == 1 and not p.is_absolute()):
        return model_path
    method = str(cfg.get("segmentation.nuclei.method", "auto"))
    if method not in ("auto", "cellpose"):
        return model_path              # the segment stage runs the classical watershed for every other value
    from .segment.backends import UnknownBackendError, cellpose_backend, get_backend

    try:
        if get_backend(cfg.get("segmentation.nuclei.backend", "cellpose")) is not cellpose_backend:
            return model_path          # classical, nuclei3d, a plugin, or a stand-in registered under the name
    except UnknownBackendError:
        return model_path              # the segment stage refuses the name itself
    if method == "auto" and not cellpose_installed():
        return model_path              # no Cellpose, so no stock weights: method auto runs the watershed, as before
    from .models_registry import MODELS, models_dir

    root = os.environ.get("GERMQUANT_MODELS")
    # an inherited relative value that load_config made absolute against its base's folder (Config.anchors)
    anchor = ((getattr(cfg, "anchors", None) or {}).get("paths") or {}).get("segmentation.nuclei.cellpose_model") or {}
    inherited = anchor.get("relative") if anchor.get("absolute") == str(model_path) else None
    rel = Path(inherited) if inherited else p
    held = Path(*rel.parts[1:]) if len(rel.parts) > 1 and rel.parts[0] == "models" else rel  # under $GERMQUANT_MODELS
    if p.is_absolute():
        looked = [p]
    else:
        looked = ([Path(root) / p] + ([Path(root) / held] if held != p else [])) if root else []
        looked += [cfg.base_dir / p, Path.cwd() / p]
    looked = list(dict.fromkeys(looked))          # the config's folder and the working directory can be one place
    folder = next((x for x in looked if x.is_dir()), None)
    msg = (f"segmentation.nuclei.cellpose_model '{model_path}' is not a file"
           + (f" ({folder} is a folder)" if folder is not None else "") + "; looked for ")
    msg += "; ".join(str(x) for x in looked)
    if inherited:
        msg += (f" (inherited as '{inherited}' from a base config in another folder and made absolute against that "
                f"folder, so GERMQUANT_MODELS is not searched)")
    elif not p.is_absolute() and not root:
        msg += " (GERMQUANT_MODELS is not set)"
    if p.is_dir():                                # Cellpose itself would be handed a folder, which it cannot load
        msg += (". Cellpose would fail to load a folder, and under method auto the run would go on with the "
                "classical watershed. Fix: ")
    else:
        msg += (". Cellpose would not stop here: it would log 'pretrained model ... not found, using default model' "
                "and segment with its stock weights. Fix: ")
    store = f"set GERMQUANT_MODELS to a folder that holds {held.as_posix()}"
    dest, n, hint = models_dir() / p.name, len(held.parts), ""    # dest: where `germquant get-model` keeps it
    download = False
    if dest.parts[-n:] == held.parts:                     # a place the lookup can be pointed at
        if dest.is_file():
            hint = f" ({Path(*dest.parts[:-n])} does)"
        elif p.name in MODELS:                            # a release asset; a local name is verified, not downloaded
            hint, download = f" (`germquant get-model {p.name}` downloads it to {dest})", True
    if not p.is_absolute() and download and dest in looked:   # the download lands where this lookup already looks
        msg += f"run `germquant get-model {p.name}` (it downloads to {dest}), or {store}"
    elif not p.is_absolute():
        msg += f"{store}{hint}, or give the model file's absolute path"
    elif inherited:
        msg += f"put the model file there, or write cellpose_model: {inherited} in this config and {store}{hint}"
    else:
        msg += f"put the model file at that path{hint}, or correct the path"
    raise FileNotFoundError(msg + ". To segment with Cellpose's stock weights on purpose, set "
                                  "segmentation.nuclei.cellpose_model: cpsam.")
