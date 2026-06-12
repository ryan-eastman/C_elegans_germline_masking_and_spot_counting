"""Resolve .nd2 channels (named by laser line) to biological roles.

Roles drive the rest of the pipeline: `dna` -> nuclei + zoning, `central_element`
-> SC tracing, `foci` -> RAD-51 counting, `axis` -> no-DAPI zoning, `granule` -> 3D objects.
"""
from __future__ import annotations

import dataclasses
from pathlib import Path

import yaml


@dataclasses.dataclass
class RoleSpec:
    names: list[str]
    index: int | None = None
    required: bool = False
    marker: str = ""


@dataclasses.dataclass
class ChannelMap:
    roles: dict[str, RoleSpec]
    match_by: str = "name"

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ChannelMap":
        data = yaml.safe_load(Path(path).read_text())
        roles = {
            role: RoleSpec(
                names=[str(n) for n in spec.get("names", [])],
                index=spec.get("index"),
                required=bool(spec.get("required", False)),
                marker=str(spec.get("marker", role)),
            )
            for role, spec in data.get("roles", {}).items()
        }
        return cls(roles=roles, match_by=data.get("match_by", "name"))

    def resolve(self, channel_names: list[str]) -> tuple[dict[str, int | None], list[str]]:
        """Return ({role: channel_index_or_None}, [qc_flags]).

        Name match is case-insensitive substring (so '477' matches 'Widefield 477').
        Falls back to the configured index when the name is missing.
        """
        lowered = [str(c).lower() for c in channel_names]
        n = len(channel_names)
        mapping: dict[str, int | None] = {}
        flags: list[str] = []

        for role, spec in self.roles.items():
            idx: int | None = None
            if self.match_by == "name":
                for want in spec.names:
                    w = want.lower()
                    hit = next((i for i, c in enumerate(lowered) if w in c or c in w), None)
                    if hit is not None:
                        idx = hit
                        break
            if idx is None and spec.index is not None and spec.index < n:
                idx = spec.index
                flags.append(f"channel:{role}:name_unmatched_used_index_{spec.index}")
            if idx is None and spec.required:
                flags.append(f"channel:{role}:MISSING_REQUIRED")
            mapping[role] = idx
        return mapping, flags

    def marker(self, role: str) -> str:
        spec = self.roles.get(role)
        return spec.marker if spec else role

    def has_role(self, role: str) -> bool:
        return role in self.roles
