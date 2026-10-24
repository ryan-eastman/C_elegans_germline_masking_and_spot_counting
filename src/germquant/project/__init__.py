"""Project folders: the manifest, recipe edits, images, models and runs of one study (docs/PROJECTS.md,
docs/DESIGN_system_packs.md, section 7).

A project folder holds germquant_project.yaml (the manifest), the pinned pack version copied into packs/, and,
only after an edit, recipe.yaml (and channel_map.yaml after a role edit). Every function takes the project
explicitly (a Project, or its folder), so a GUI calling from a worker thread never depends on a context it
did not set.

Importing this package loads nothing heavy (no torch, CuPy, SpotMAX, Cellpose, bioio or dask.array): it is
the module the torch-free entry ``python -m germquant.project`` runs. Part 1a of unit P5-PRJ: the frozen
signatures of PLAN.md section 9.1, as stubs that raise NotImplementedError naming the unit that builds them.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "Project", "add_images", "break_lock", "collect_project", "export_commands", "list_rounds", "new_project",
    "new_round", "open_project", "preview_roles", "project_config", "project_status", "revert_recipe",
    "run_project", "save_project", "use_model", "write_recipe",
]


def _todo(name: str, unit: str) -> str:
    return (f"germquant.project.{name} is not built yet: unit {unit} of tranche 5 Train B builds it "
            f"(C:/Users/ryane/germquant_golden/t5_train_b/PLAN.md, section 9.1)")


@dataclass
class Project:
    """An open project: ``root`` (the folder), ``manifest`` (germquant_project.yaml as a dict), ``warnings``
    (what open_project found), and the manifest's ``images`` and ``pack`` entries as attributes."""

    root: Path
    manifest: dict
    warnings: list[str] = field(default_factory=list)
    config_hash: str | None = None

    @property
    def images(self) -> list[dict]:
        return self.manifest.setdefault("images", [])

    @property
    def pack(self) -> dict:
        return self.manifest.setdefault("pack", {})


def new_project(root, pack, recipe=None, images=(), *, name=None, plugin_version=None) -> Project:
    """Create a project in ``root``: copy the pinned pack version (``pack`` is ``id`` or ``id@version``) into
    <root>/packs/, write the manifest and no recipe.yaml, and add ``images``."""
    raise NotImplementedError(_todo("new_project", "P5-PRJ"))


def open_project(root) -> Project:
    """Read the manifest of the project in ``root`` and recompute its config hash (a mismatch is a warning)."""
    raise NotImplementedError(_todo("open_project", "P5-PRJ"))


def save_project(project, *, plugin_version=None, device=None, interpreters=None, runner=None, xy_stride=None,
                 z_range=None, overrides=None, samples=None) -> Project:
    """Write the run settings and the other non-recipe fields of the manifest (None leaves a field as it is)."""
    raise NotImplementedError(_todo("save_project", "P5-PRJ"))


def add_images(project, paths, *, plugin_version=None) -> dict:
    """Add files, folders and globs: ``{"added": [image_id], "refused": [{"path": str, "reason": str}]}``."""
    raise NotImplementedError(_todo("add_images", "P5-PRJ"))


def project_config(project):
    """The Config the project runs with (cfg.pack_record set from the manifest, with "edited")."""
    raise NotImplementedError(_todo("project_config", "P5-PRJ"))


def preview_roles(project, roles, *, image_id=None) -> dict:
    """``{"yaml": str, "bindings": {image_id: {role: {"channel", "how"}}}, "problems": [str]}``: what a role edit
    would write, without writing it."""
    raise NotImplementedError(_todo("preview_roles", "P5-PRJ"))


def write_recipe(project, *, edits=None, unset=(), roles=None, origin="form", plugin_version=None) -> str:
    """Write recipe.yaml (and channel_map.yaml for role edits), validated by origin; returns the config_hash."""
    raise NotImplementedError(_todo("write_recipe", "P5-PRJ"))


def revert_recipe(project, *, plugin_version=None) -> str:
    """Remove every recipe edit (no recipe.yaml): the project hashes as the pack recipe again."""
    raise NotImplementedError(_todo("revert_recipe", "P5-PRJ"))


def use_model(project, model, *, card=None, plugin_version=None) -> str:
    """Name a model of the project's training/ folder by relative path in recipe.yaml; returns the config_hash."""
    raise NotImplementedError(_todo("use_model", "P5-PRJ"))


def run_project(project, images=None, resume=True, device=None, shard=None, progress=None, *, via="cli",
                run_id=None, should_stop=None) -> dict:
    """Run the project's images: ``{"run_id", "log", "counts": {"done", "failed", "skipped"}, "exit_code"}``."""
    raise NotImplementedError(_todo("run_project", "P5-RUN"))


def project_status(project) -> dict:
    """The status record of PLAN.md section 9.3 (what `germquant project status --json` prints)."""
    raise NotImplementedError(_todo("project_status", "P5-RUN"))


def collect_project(project) -> dict:
    """Collect the per-image tables of results/ with the recipe's qc.germline_outlier_factor."""
    raise NotImplementedError(_todo("collect_project", "P5-RUN"))


def export_commands(project, *, shell="sh", shards=None) -> str:
    """The commands that run the project headless (sh or ps1), optionally split into shards."""
    raise NotImplementedError(_todo("export_commands", "P5-RUN"))


def break_lock(project, *, host, pid) -> bool:
    """Remove the project's lock file when its holder is ``host`` and ``pid``."""
    raise NotImplementedError(_todo("break_lock", "P5-PRJ"))


def list_rounds(project) -> list[str]:
    """The training rounds under training/."""
    raise NotImplementedError(_todo("list_rounds", "P5-PRJ"))


def new_round(project, name=None, *, plugin_version=None) -> Path:
    """Create training/<name>/ (default the next of r01, r02 and so on) and return it."""
    raise NotImplementedError(_todo("new_round", "P5-PRJ"))
