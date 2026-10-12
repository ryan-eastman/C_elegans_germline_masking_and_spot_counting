"""Helpers for pack maintainers: the frozen base of an outside pack (docs/DESIGN_system_packs.md, section 3.12).

Every outside recipe extends its pack's own ``config/profiles/pack_base.yaml`` instead of generic.yaml (which
extends config.yaml), so the planned config.yaml hash bump never moves an outside pack's hash. The base is
written ONCE from a resolved profile (yaml.safe_dump, ASCII, no comments), with io.channel_map pointing at the
pack's own map, no machine path, every worm stage switch off, the neutral filename regex and the QC settings
of a field that is not a gonad, and without the keys on the lint deny-list (spots.channel and the other keys
nothing reads). After that it is frozen text: lint and the recipe goldens catch any edit, and write_frozen_base
refuses to replace an existing base with other text unless told to rebase (which moves every recipe hash of the
pack, so it goes with a minor version bump and a CHANGELOG entry).
"""
from __future__ import annotations

import copy
import difflib
from pathlib import Path

import yaml

from .lint import DENY_KEYS, NEUTRAL_REGEX, WORM_SWITCHES


class FrozenBaseError(ValueError):
    """An existing frozen pack_base.yaml that a write would replace with other text (the message holds the diff)."""


def _set(d: dict, dotted: str, value) -> None:
    parts = dotted.split(".")
    for p in parts[:-1]:
        d = d.setdefault(p, {})
    d[parts[-1]] = value


def _drop(d: dict, dotted: str) -> None:
    parts = dotted.split(".")
    for p in parts[:-1]:
        d = d.get(p)
        if not isinstance(d, dict):
            return
    d.pop(parts[-1], None)


def frozen_base_data(profile: str | Path, *, channel_map: str, germline_outlier_factor: float = 1.0e9,
                     min_nuclei: int = 10, overrides: dict | None = None) -> dict:
    """The data of a frozen pack base: ``profile`` resolved through load_config, with the pack's channel map,
    every worm switch off, the neutral regex, the outside QC settings and the deny-listed keys dropped;
    ``overrides`` (dotted keys) are applied last."""
    from ..config import load_config

    data = copy.deepcopy(load_config(Path(profile)).as_dict())
    for key in DENY_KEYS:
        _drop(data, key)
    _set(data, "io.channel_map", str(channel_map))
    for sw in WORM_SWITCHES:
        _set(data, sw, False)
    _set(data, "metadata.filename_regex", NEUTRAL_REGEX)
    _set(data, "qc.germline_outlier_factor", float(germline_outlier_factor))
    _set(data, "qc.min_nuclei", int(min_nuclei))
    for key, value in (overrides or {}).items():
        _set(data, key, value)
    return data


def frozen_base_text(data: dict) -> str:
    """The base as frozen text: yaml.safe_dump, keys sorted, ASCII only, LF."""
    text = yaml.safe_dump(data, sort_keys=True, default_flow_style=False, allow_unicode=False)
    text.encode("ascii")                               # raises if anything slipped through
    return text


def base_diff(path: str | Path, text: str, limit: int = 80) -> str:
    """A unified diff from the base at ``path`` (LF-normalised) to ``text``, at most ``limit`` lines."""
    old = Path(path).read_bytes().replace(b"\r\n", b"\n").decode("ascii", "replace")
    lines = list(difflib.unified_diff(old.splitlines(), text.splitlines(), str(path), "regenerated", lineterm=""))
    return "\n".join(lines[:limit] + ([f"... {len(lines) - limit} more lines"] if len(lines) > limit else []))


def base_unchanged(path: str | Path, data: dict) -> bool:
    """Whether the base at ``path`` already holds exactly the text ``data`` makes (LF-normalised)."""
    p = Path(path)
    return p.is_file() and p.read_bytes().replace(b"\r\n", b"\n") == frozen_base_text(data).encode("ascii")


def write_frozen_base(path: str | Path, data: dict, *, rebase: bool = False) -> Path:
    """Write the frozen base to ``path`` and return it. An existing base with the same text (LF-normalised) is
    left as it is. An existing base with other text raises FrozenBaseError with the diff, unless ``rebase``:
    that moves every recipe hash of the pack, so it needs a minor version bump and a CHANGELOG entry (design 3.7
    and 3.12)."""
    p = Path(path)
    text = frozen_base_text(data)
    if p.is_file():
        if base_unchanged(p, data):
            return p
        if not rebase:
            raise FrozenBaseError(f"{p} is frozen and the regenerated text differs; writing it would move every "
                                  f"recipe hash of the pack. Pass --rebase only together with a minor version bump "
                                  f"and a CHANGELOG entry.\n{base_diff(p, text)}")
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(text.encode("ascii"))
    return p
