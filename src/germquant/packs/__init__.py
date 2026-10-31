"""System packs: checksummed folders of static files that a user picks by name (docs/PACKS.md,
docs/DESIGN_system_packs.md, section 3).

Choosing a pack and one of its recipes resolves to a profile file that load_config reads like any other, so
``--pack X --recipe R`` and ``--config <that file>`` give the same config_hash and the same numbers. A pack
names backends and models by name, never by callable; weights live once per digest in the model store.

Packs are data. Importing this package, listing packs and the shallow pack_status parse YAML and JSON with
yaml.safe_load and never import a backend, torch, Cellpose, SpotMAX or CuPy (gate G6), so a GUI can refresh a
pack card cheaply and the cuBLAS load order on the RTX 5090 (torch before CuPy) is never disturbed. Code that
must stay light imports this package, not germquant.api (which imports the pipeline and with it torch). The
torch-free command entry is ``python -m germquant.packs`` (list, show, check, fetch-models, fetch-example,
install, verify, build-sums, lint, derive, --version; ``--project DIR`` on list, show, check and fetch-models
searches the project's packs folder first). No code from a pack folder ever runs.
"""
from __future__ import annotations

import re
from pathlib import Path

from ..models_registry import ModelPinError, check_model_pin, resolve_model_name
from .blobs import blobs_dir, check_blob, fetch_models, find_model_spec, model_location, resolve_pack_model
from .derive import PackDeriveError, derive_pack, portable_card
from .discovery import (
    Discovery,
    discover,
    engine_compatible,
    engine_root_for,
    env_pack_dirs,
    find_pack,
    install_root,
    pack_roots,
    parse_ref,
    recipe_profile_path,
    refresh,
    resolve_recipe,
    use_project,
    user_packs_dir,
)
from .examples import ExampleError, check_expect, examples_dir, fetch_example, find_example, list_examples
from .format import (
    MANIFEST,
    SUMS_FILE,
    ModelSpec,
    Pack,
    PackFormatError,
    Recipe,
    build_sums_text,
    load_pack_folder,
    sums_sha256,
    validate_manifest,
    verify_sums,
    write_sums,
)
from .install import PackInstallError, install_pack, vendor
from .lint import lint_pack
from .status import engine_files_problems, pack_status
from .validate import tunable_problems, validate_pack

__all__ = [
    "MANIFEST", "SUMS_FILE", "Discovery", "ExampleError", "ModelPinError", "ModelSpec", "Pack", "PackFormatError",
    "PackInstallError", "Recipe", "check_model_pin", "resolve_model_name",
    "blobs_dir", "build_sums", "build_sums_text", "check_blob", "check_expect", "discover", "engine_compatible",
    "engine_files_problems", "engine_root_for", "env_pack_dirs", "examples_dir", "fetch_example", "fetch_models",
    "fetch_pack_models", "find_example", "find_model_spec", "find_pack", "install_pack", "install_root",
    "lint_pack", "list_examples", "list_packs", "load_pack",
    "load_pack_folder", "model_location", "no_builtin_root_message", "pack_record", "pack_roots", "pack_status",
    "parse_ref", "recipe_profile_path", "refresh", "resolve_pack_model", "resolve_recipe", "sums_sha256",
    "tunable_problems", "use_project", "user_packs_dir", "validate_manifest", "validate_pack", "vendor",
    "verify_pack", "verify_sums", "version_info", "version_line", "write_sums",
    "FEATURE_PROBES", "PackDeriveError", "derive_pack", "features", "portable_card",
]


def list_packs(all: bool = False, project: str | Path | None = None) -> list[dict]:  # noqa: A002 - the CLI flag
    """The packs discovery finds, as dicts (Pack.summary plus ``state``): the usable ones in search order, and
    with ``all`` also the unusable ones (with their problems) and the shadowed ones."""
    disc = discover(project)
    out = [{**p.summary(), "state": "ok"} for p in disc.packs]
    if all:
        out += [{"id": u["id"], "version": u["version"], "folder": u["folder"], "source": u["source"],
                 "problems": u["problems"], "state": "unusable", "recipes": [], "models": []} for u in disc.unusable]
        out += [{**p.summary(), "state": "shadowed"} for p in disc.shadowed]
    return out


def no_builtin_root_message() -> str | None:
    """The message `pack list` prints when no built-in root exists (a non-editable install), or None."""
    from ..config import builtin_roots

    if builtin_roots():
        return None
    return ("No built-in packs: this germquant is not an editable checkout (no src/germquant above the installed "
            "package) and GERMQUANT_CONFIG_ROOT is not set. Two fixes: install the engine as an editable checkout "
            "(from the checkout: uv pip install -e .), the supported route for 0.2.0, or set GERMQUANT_CONFIG_ROOT "
            "to a checkout that holds config/ and packs/, as the container does.")


def load_pack(ref_or_folder) -> Pack:
    """A Pack from a folder holding pack.yaml, or from a name (``id`` or ``id@version``)."""
    if isinstance(ref_or_folder, Pack):
        return ref_or_folder
    p = Path(str(ref_or_folder))
    if (p / MANIFEST).is_file():
        return load_pack_folder(p)
    return find_pack(str(ref_or_folder))


def pack_record(pack, recipe: str | None = None) -> dict:
    """What a pack run records in run_geometry["pack"] (and the batch-root manifest): id, version, recipe,
    source, pack_sha256, object_kind, dims and status. Never a table column, the done marker or a hash."""
    p = load_pack(pack)
    r = p.recipe(recipe)
    return {"id": p.id, "version": p.version, "recipe": r.id, "source": p.source,
            "pack_sha256": sums_sha256(p.folder), "object_kind": p.object_kind, "dims": r.dims, "status": r.status}


def fetch_pack_models(pack, recipe: str | None = None, from_files: dict | None = None, *, download: bool = True,
                      say=None) -> list[dict]:
    """`germquant pack fetch-models`: put the models a recipe (or the whole pack) needs into the store."""
    return fetch_models(load_pack(pack), recipe, from_files, download=download, say=say)


def build_sums(folder: str | Path) -> Path:
    """`germquant pack build-sums DIR`: (re)write DIR/SHA256SUMS."""
    return write_sums(Path(folder))


def verify_pack(ref_or_folder) -> list[str]:
    """`germquant pack verify`: problems with the manifest, SHA256SUMS and (root: checkout) the engine files."""
    try:
        p = load_pack(ref_or_folder)
    except (PackFormatError, FileNotFoundError) as e:
        return [str(e)]
    return [f"manifest: {x}" for x in p.problems] + verify_sums(p.folder) + engine_files_problems(p)


# The features `python -m germquant.packs --version --json` lists (docs/PACKS.md, "Versions and features"), in
# the order of the Train B contract. The version string stays 0.2.0.dev0 through tranche 5, so a client such as
# the napari plugin gates on these names. A feature is listed when its code is in this engine: its module file
# exists and, where a name is given (for planar, projects, train_curve and post_stages, whose module or package
# exists before the feature does), the module defines that name at top level. projects waits for project/run.py and
# its run_project: the package's first stubs (every call raising NotImplementedError) are not the feature; so does
# post_stages wait for the run commands to take post-stages, marked by api.py's list_post_stages (the
# germquant.post_stages package of phase CORE alone is not the feature). Files are read, never imported (gate G6),
# so an engine that lacks a feature never claims it.
FEATURE_PROBES = (
    ("packs", "packs/__init__.py", None),
    ("planar", "planar.py", "PLANAR_DZ"),
    ("labels_from", "io/labels.py", None),
    ("projects", "project/run.py", "run_project"),
    ("sample_sheet", "io/sample_sheet.py", None),
    ("import_labels", "import_labels.py", None),
    ("train_curve", "segment/training.py", "train_curve"),
    ("pack_derive", "packs/derive.py", None),
    ("tune_spots", "tune/__init__.py", None),
    ("post_stages", "api.py", "list_post_stages"),
)


def _defines(path: Path, name: str) -> bool:
    """Whether the Python file ``path`` defines ``name`` at top level (a def, a class or an assignment)."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    n = re.escape(name)
    return re.search(rf"^(?:def|class)\s+{n}\b|^{n}\s*(?::[^=\n]*)?=(?!=)", text, re.M) is not None


def features(package_dir: str | Path | None = None) -> list[str]:
    """The features of this engine (FEATURE_PROBES), or those a germquant package folder ``package_dir`` holds."""
    root = Path(package_dir) if package_dir is not None else Path(__file__).resolve().parents[1]
    return [name for name, rel, attr in FEATURE_PROBES
            if (root / rel).is_file() and (attr is None or _defines(root / rel, attr))]


def version_info() -> dict:
    """The engine version, git sha, python and interpreter (what `--version` prints; the plugin compares it),
    and the features of this engine (FEATURE_PROBES)."""
    import platform
    import sys

    from .. import __version__
    from ..provenance import git_sha

    return {"germquant": __version__, "git_sha": git_sha(), "python": platform.python_version(),
            "executable": sys.executable, "features": features()}


def version_line() -> str:
    v = version_info()
    return f"germquant {v['germquant']} git {v['git_sha']} python {v['python']}"
