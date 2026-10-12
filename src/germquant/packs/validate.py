"""packs.validate_pack: what the loader refuses (docs/DESIGN_system_packs.md, section 3.3).

validate_manifest (packs.format) checks pack.yaml on its own: the format major, the id, recipe and model name
patterns, model digests and pinned names, the recipe layout and validated claims. validate_pack adds what
needs the recipes loaded: a tunable key that the recipe does not set, or whose value is not a scalar. The
parameter form is built from the tunables, so such a knob would add a key, and with it a new config hash,
instead of changing a value the recipe already sets.

Who uses it: install refuses a pack that validate_pack faults (packs.install.check_installable); pack status
blocks a recipe whose own tunables are faulty; lint reports every fault. Discovery stays manifest-only, so
listing packs never loads a recipe. A recipe that does not load is not reported here (pack status names it as
a blocker and lint as a problem), so installing packs that extend each other does not depend on their order.
"""
from __future__ import annotations

from pathlib import Path

from .format import Pack, Recipe, load_pack_folder

_MISSING = object()


def _lookup(d: dict, dotted: str):
    node = d
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return _MISSING
        node = node[part]
    return node


def tunable_problems(recipe: Recipe, merged: dict) -> list[str]:
    """Problems with a recipe's tunables against its resolved config (``Config.as_dict()``): each key must be
    set by the recipe chain and hold a scalar (a mapping or a list is not one)."""
    out = []
    for t in recipe.tunables or []:
        key = str(t.get("key", "") or "") if isinstance(t, dict) else ""
        if not key:
            continue                                   # a manifest problem (validate_manifest names it)
        value = _lookup(merged, key)
        if value is _MISSING:
            out.append(f"recipe {recipe.id}: tunable {key} is not set by the recipe")
        elif isinstance(value, (dict, list)):
            out.append(f"recipe {recipe.id}: tunable {key} is not a scalar")
    return out


def load_recipe_config(pack: Pack, recipe):
    """The Config of a recipe (a Recipe or a recipe id), read through load_config like any profile."""
    from ..config import load_config
    from .discovery import recipe_profile_path

    return load_config(recipe_profile_path(pack, recipe))


def validate_pack(pack_or_folder) -> list[str]:
    """Every problem the loader refuses a pack for (module docstring): validate_manifest's, then the tunables of
    each recipe that loads. ``pack_or_folder`` is a Pack or a folder holding pack.yaml. Empty when the pack is
    sound; PackFormatError when pack.yaml cannot be parsed at all. No recipe is loaded while the manifest has
    problems (a profile path could point outside the pack)."""
    pack = pack_or_folder if isinstance(pack_or_folder, Pack) else load_pack_folder(Path(pack_or_folder))
    probs = list(pack.problems)
    if probs:
        return probs
    for r in pack.recipes:
        if not r.tunables:
            continue
        try:
            cfg = load_recipe_config(pack, r)
        except Exception:  # noqa: BLE001 - an unloadable recipe is pack status's and lint's to report
            continue
        probs += tunable_problems(r, cfg.as_dict())
    return probs
