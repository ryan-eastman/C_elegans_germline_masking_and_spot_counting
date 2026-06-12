"""Config loading. Thin wrapper over YAML with attribute access + the channel map."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import yaml

from .io.channel_map import ChannelMap


class Config:
    """Dict-backed config with attribute access and a resolved ChannelMap.

    Access nested keys via attributes (``cfg.segmentation.nuclei.method``) or
    ``cfg.get("segmentation.nuclei.method", default)``.
    """

    def __init__(self, data: dict, *, base_dir: Path, channel_map: ChannelMap, raw_text: str):
        self._data = data
        self.base_dir = base_dir
        self.channel_map = channel_map
        self._raw_text = raw_text

    # ---- access helpers ----
    def __getattr__(self, name: str) -> Any:
        try:
            val = self._data[name]
        except KeyError as e:  # pragma: no cover
            raise AttributeError(name) from e
        return _Node(val) if isinstance(val, dict) else val

    def get(self, dotted: str, default: Any = None) -> Any:
        node: Any = self._data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    @property
    def hash(self) -> str:
        """Stable hash of the resolved config — goes into the provenance manifest."""
        return hashlib.sha256(self._raw_text.encode()).hexdigest()[:12]

    def as_dict(self) -> dict:
        return json.loads(json.dumps(self._data))


class _Node:
    def __init__(self, d: dict):
        self._d = d

    def __getattr__(self, name: str) -> Any:
        try:
            val = self._d[name]
        except KeyError as e:
            raise AttributeError(name) from e
        return _Node(val) if isinstance(val, dict) else val

    def get(self, key: str, default: Any = None) -> Any:
        return self._d.get(key, default)

    def __repr__(self) -> str:  # pragma: no cover
        return f"_Node({self._d!r})"


def load_config(path: str | Path) -> Config:
    path = Path(path)
    raw_text = path.read_text()
    data = yaml.safe_load(raw_text)
    base_dir = path.parent.parent if path.parent.name == "config" else path.parent

    cm_rel = data.get("io", {}).get("channel_map")
    if cm_rel is None:
        raise ValueError("config.io.channel_map is required")
    cm_path = (base_dir / cm_rel) if not Path(cm_rel).is_absolute() else Path(cm_rel)
    channel_map = ChannelMap.from_yaml(cm_path)

    return Config(data, base_dir=base_dir, channel_map=channel_map, raw_text=raw_text)
