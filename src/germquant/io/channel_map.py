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

        Name match (match_by=='name') is a case-insensitive substring test: the configured
        token must appear *in* the actual channel name (e.g. '477' matches 'Widefield 477').
        The match is one-directional on purpose — an earlier version also tested ``c in w``,
        which made an empty/blank channel name ('' in 'dapi' == True) silently match the first
        role and swap DAPI/SYP/RAD-51. Blank channel names are skipped here for the same reason.
        Falls back to the configured ``index`` when the name doesn't match (or in index mode).
        """
        lowered = [str(c).lower() for c in channel_names]
        n = len(channel_names)
        mapping: dict[str, int | None] = {}
        flags: list[str] = []

        for role, spec in self.roles.items():
            idx: int | None = None
            if self.match_by == "name":
                for want in spec.names:
                    w = want.lower().strip()
                    if not w:
                        continue
                    hit = next((i for i, c in enumerate(lowered) if c and w in c), None)
                    if hit is not None:
                        idx = hit
                        break
                if idx is None and spec.index is not None and spec.index < n:
                    idx = spec.index
                    flags.append(f"channel:{role}:name_unmatched_used_index_{spec.index}")
            else:  # match_by == "index": the configured index is the primary resolution
                if spec.index is not None and spec.index < n:
                    idx = spec.index
            if idx is None and spec.required:
                flags.append(f"channel:{role}:MISSING_REQUIRED")
            mapping[role] = idx

        # two roles must never resolve to the same channel (they'd read identical pixels)
        seen: dict[int, str] = {}
        for role, i in mapping.items():
            if i is None:
                continue
            if i in seen:
                flags.append(f"channel:role_collision:{seen[i]}={role}=idx{i}")
            else:
                seen[i] = role
        return mapping, flags

    def marker(self, role: str) -> str:
        spec = self.roles.get(role)
        return spec.marker if spec else role

    def has_role(self, role: str) -> bool:
        return role in self.roles
