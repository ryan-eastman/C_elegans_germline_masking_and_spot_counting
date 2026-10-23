"""Imported SC geometry: Imaris Filaments (.ims, read with h5py), SWC and a CSV of paths.

Every path comes back as unbranched point lists in the pipeline frame (um, z y x). Imaris positions are voxel
centres in stage microns: ``index = (p - ExtMin) / v_imaris - 0.5`` (continuous), then ``um = index x our voxel``
(the run's voxel; for a TIFF without a z spacing it is given explicitly). A branched filament is split at its
branch points; single-vertex filaments are ignored. The stored Imaris ``Segment Length`` statistics are kept for
the as-drawn check of the validation protocol (S3a).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ...fsutil import long_path


@dataclass
class ImportedPath:
    path_id: str
    points: np.ndarray              # (n, 3) um, pipeline frame (z, y, x)
    verified: bool = False
    time_index: int | None = None
    source: str = ""
    stored_length_um: float | None = None
    extra: dict = field(default_factory=dict)


def _paths_from_edges(n: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    """Unbranched vertex sequences of a graph (split at nodes whose degree is not 2; cycles kept whole)."""
    adj: dict[int, list[int]] = {i: [] for i in range(n)}
    for a, b in edges:
        if a == b:
            continue
        adj[a].append(b)
        adj[b].append(a)
    used: set[tuple[int, int]] = set()
    out = []

    def walk(start, nxt):
        seq = [start, nxt]
        used.add((min(start, nxt), max(start, nxt)))
        prev, cur = start, nxt
        while len(adj[cur]) == 2:
            nb = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            e = (min(cur, nb), max(cur, nb))
            if e in used:
                break
            used.add(e)
            seq.append(nb)
            prev, cur = cur, nb
        return seq

    for v in sorted(adj):
        if len(adj[v]) != 2:
            for nb in sorted(adj[v]):
                e = (min(v, nb), max(v, nb))
                if e not in used:
                    out.append(walk(v, nb))
    for v in sorted(adj):                       # pure cycles
        for nb in sorted(adj[v]):
            e = (min(v, nb), max(v, nb))
            if e not in used:
                out.append(walk(v, nb))
    return [p for p in out if len(p) >= 2]


def imaris_to_um(pos_xyz: np.ndarray, ext_min_xyz, v_imaris_xyz, our_voxel_zyx) -> np.ndarray:
    idx = (np.asarray(pos_xyz, float) - np.asarray(ext_min_xyz, float)) / np.asarray(v_imaris_xyz, float) - 0.5
    um_xyz = idx * np.asarray(our_voxel_zyx, float)[::-1]
    return um_xyz[:, ::-1]


def read_ims_filaments(path: str | Path, our_voxel_zyx, time_index: int | None = None) -> list[ImportedPath]:
    import h5py

    from ...validate.imaris_ims import _content_group, _decode, read_geometry

    out: list[ImportedPath] = []
    with h5py.File(long_path(path), "r") as f:
        geom = read_geometry(f)
        content = _content_group(f)
        if content is None:
            return out
        names = sorted((k for k in content.keys() if k.startswith("Filaments") and k[9:].isdigit()),
                       key=lambda k: int(k[9:]))
        for gname in names:
            g = content[gname]
            if "Vertex" not in g or "Filament" not in g:
                continue
            V = g["Vertex"][()]
            E = g["Edge"][()] if "Edge" in g else np.zeros(0, dtype=[("VertexA", "<i8"), ("VertexB", "<i8")])
            Fl = g["Filament"][()]
            stats = {}
            if "StatisticsType" in g and "StatisticsValue" in g:
                st = g["StatisticsType"][()]
                names_by_id = {int(r["ID"]): _decode(r["Name"]) for r in st}
                for r in g["StatisticsValue"][()]:
                    nm = names_by_id.get(int(r["ID_StatisticsType"]), "")
                    stats.setdefault(nm, {})[int(r["ID_Object"])] = float(r["Value"])
            pos = np.stack([V["PositionX"], V["PositionY"], V["PositionZ"]], 1).astype(float)
            for fl in Fl:
                t = int(fl["ID_Time"]) if "ID_Time" in Fl.dtype.names else 0
                if time_index is not None and t != int(time_index):
                    continue
                vb, ve = int(fl["IndexVertexBegin"]), int(fl["IndexVertexEnd"])
                eb, ee = int(fl["IndexEdgeBegin"]), int(fl["IndexEdgeEnd"])
                nv = ve - vb
                if nv < 2:
                    continue                       # single-vertex stubs
                ed = [(int(a), int(b)) for a, b in zip(E["VertexA"][eb:ee], E["VertexB"][eb:ee])]
                if ed and all(vb <= a < ve and vb <= b < ve for a, b in ed):
                    ed = [(a - vb, b - vb) for a, b in ed]            # global indices
                elif ed and not all(0 <= a < nv and 0 <= b < nv for a, b in ed):
                    raise ValueError(f"{path}: {gname} filament {int(fl['ID'])}: edge indices outside its vertices")
                um = imaris_to_um(pos[vb:ve], geom["ext_min_xyz"], geom["voxel_xyz_um"], our_voxel_zyx)
                for k, seq in enumerate(_paths_from_edges(nv, ed)):
                    out.append(ImportedPath(path_id=f"{gname}:{int(fl['ID'])}:{k}", points=um[seq], time_index=t,
                                            source=str(path),
                                            stored_length_um=stats.get("Filament Length (sum)", {}).get(int(fl["ID"])),
                                            extra={"n_vertices": nv}))
    return out


def read_swc(path: str | Path, spacing_zyx=None) -> list[ImportedPath]:
    """SWC (id type x y z radius parent); coordinates in um, or in voxels when ``spacing_zyx`` is given."""
    rows = []
    with open(long_path(path), encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            rows.append((int(parts[0]), float(parts[2]), float(parts[3]), float(parts[4]), int(parts[6])))
    if not rows:
        return []
    ids = [r[0] for r in rows]
    index = {i: k for k, i in enumerate(ids)}
    xyz = np.array([[r[1], r[2], r[3]] for r in rows])
    zyx = xyz[:, ::-1] * (np.asarray(spacing_zyx, float) if spacing_zyx is not None else 1.0)
    edges = [(index[r[0]], index[r[4]]) for r in rows if r[4] in index]
    return [ImportedPath(path_id=f"swc:{k}", points=zyx[seq], source=str(path))
            for k, seq in enumerate(_paths_from_edges(len(rows), edges))]


def read_paths_csv(path: str | Path, image_id: str | None = None) -> list[ImportedPath]:
    """Columns image_id (optional), strand_id, point_index, z_um, y_um, x_um, verified (optional)."""
    df = pd.read_csv(long_path(path))
    if image_id is not None and "image_id" in df.columns:
        df = df[df["image_id"].astype(str) == str(image_id)]
    out = []
    for sid, g in df.groupby("strand_id", sort=True):
        g = g.sort_values("point_index")
        ver = bool(g["verified"].astype(str).str.lower().isin(["true", "1"]).all()) if "verified" in g else False
        P = g[["z_um", "y_um", "x_um"]].to_numpy(float)
        if len(P) >= 2:
            out.append(ImportedPath(path_id=f"csv:{sid}", points=P, verified=ver, source=str(path)))
    return out


def read_spec(kind: str, path: str, *, image_id: str, spacing) -> list[ImportedPath]:
    kind = kind.lower()
    if kind == "ims":
        return read_ims_filaments(path, spacing)
    if kind == "swc":
        return read_swc(path)
    if kind == "csv":
        return read_paths_csv(path, image_id)
    raise ValueError(f"--paths-from wants ims:FILE, swc:FILE or csv:FILE, got {kind}:")


def assign_nucleus(points: np.ndarray, labels: np.ndarray, spacing, search_vox=(2, 4, 4)) -> int | None:
    """The label holding most of the path's points; else the most common label near the points."""
    idx = np.rint(np.asarray(points, float) / np.asarray(spacing, float)).astype(int)
    shp = np.asarray(labels.shape)
    ok = np.all((idx >= 0) & (idx < shp), axis=1)
    if not ok.any():
        return None
    vals = labels[tuple(idx[ok].T)].astype(np.int64)
    nz = vals[vals > 0]
    if len(nz):
        return int(np.bincount(nz).argmax())
    found = []
    for p in idx[ok]:
        sl = tuple(slice(max(0, c - r), min(n, c + r + 1)) for c, r, n in zip(p, search_vox, shp))
        sub = labels[sl]
        found.extend(sub[sub > 0].tolist())
    return int(np.bincount(np.asarray(found)).argmax()) if found else None
