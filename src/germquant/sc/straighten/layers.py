"""napari layer data for sc-straighten outputs, without importing napari (LayerData tuples:
``(data, kwargs, layer_type)``). World coordinates are the pipeline frame in microns (z, y, x), so these layers
overlay the run's label image loaded with ``scale = spacing``.

* ``strand_layers``: the strands as a Shapes path layer (features = the strands table, coloured by class), the
  end points (by end class), the gap midpoints and the foci (by assignment).
* ``strip_layers``: a gallery of one nucleus's straightened strands (flat strips stacked as rows, NaN padded,
  scale (du, ds) um) and its foci at (row, u, s).
* ``image_layers``: a crop of a channel with ``scale = spacing`` and ``translate = origin x spacing``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ...fsutil import long_path
from . import schema as S
from .outputs import CLASS_COLOURS, read_strips

ASSIGN_COLOURS = {"unique": "#ffd700", "ambiguous": "#ff8c00", "off_strand": "#808080", "no_strand": "#404040"}
END_COLOURS = {"free_fit": "#2ca02c", "free_unfit": "#ff7f0e", "contact": "#1f77b4", "stack_edge": "#d62728",
               "territory_edge": "#9467bd", "cycle_cut": "#8c564b", "max_length": "#e377c2", "as_drawn": "#7f7f7f"}


def _read(results_dir, image_id, key) -> pd.DataFrame:
    p = Path(results_dir) / f"{image_id}{S.SUFFIXES[key]}"
    if not p.is_file():
        return pd.DataFrame(columns=S.TABLES.get(key, []))
    return pd.read_csv(long_path(p), low_memory=False)


def strand_layers(results_dir, image_id: str, nucleus_ids=None) -> list[tuple]:
    strands = _read(results_dir, image_id, "strands")
    samples = _read(results_dir, image_id, "samples")
    gaps = _read(results_dir, image_id, "gaps")
    foci = _read(results_dir, image_id, "foci")
    if nucleus_ids is not None:
        keep = {int(v) for v in nucleus_ids}
        strands = strands[strands["nucleus_id"].isin(keep)]
        samples = samples[samples["nucleus_id"].isin(keep)]
        gaps = gaps[gaps["nucleus_id"].isin(keep)] if len(gaps) else gaps
        foci = foci[foci["nucleus_id"].isin(keep)] if len(foci) else foci
    layers = []
    paths, colours, feats = [], [], []
    st = strands.set_index("strand_id") if len(strands) else pd.DataFrame()
    if len(samples):
        smp = samples[samples["in_strand"].astype(str).str.lower().isin(["true", "1"])]
        for sid, g in smp.groupby("strand_id", sort=True):
            g = g.sort_values("s_um")
            if len(g) < 2 or sid not in st.index:
                continue
            paths.append(g[["z_um", "y_um", "x_um"]].to_numpy(float))
            colours.append(CLASS_COLOURS.get(st.loc[sid, "strand_class"], "#ff7f0e"))
            feats.append(st.loc[sid].to_dict() | {"strand_id": int(sid)})
    feat_df = pd.DataFrame(feats)
    layers.append((paths, {"name": f"{image_id} strands", "shape_type": "path", "edge_color": colours,
                           "edge_width": 0.08, "features": feat_df, "opacity": 0.9}, "shapes"))
    if len(strands):
        ends = np.vstack([strands[["end0_z_um", "end0_y_um", "end0_x_um"]].to_numpy(float),
                          strands[["end1_z_um", "end1_y_um", "end1_x_um"]].to_numpy(float)])
        cls = list(strands["end0_class"]) + list(strands["end1_class"])
        layers.append((ends, {"name": f"{image_id} strand ends", "size": 0.15,
                              "face_color": [END_COLOURS.get(c, "#ffffff") for c in cls],
                              "features": pd.DataFrame({"end_class": cls})}, "points"))
    if len(gaps):
        layers.append((gaps[["mid_z_um", "mid_y_um", "mid_x_um"]].to_numpy(float),
                       {"name": f"{image_id} gaps", "size": 0.2, "symbol": "x", "face_color": "#ff00ff",
                        "features": gaps.reset_index(drop=True)}, "points"))
    if len(foci):
        layers.append((foci[["z_um", "y_um", "x_um"]].to_numpy(float),
                       {"name": f"{image_id} foci", "size": 0.25,
                        "face_color": [ASSIGN_COLOURS.get(a, "#ffffff") for a in foci["assignment"]],
                        "features": foci.reset_index(drop=True)}, "points"))
    return layers


def strip_gallery(results_dir, image_id: str, nucleus_id: int, kind: str = "flat", role: str = "central_element"):
    """(gallery (rows, n_u, max_len) float32 with NaN padding, strands table rows, description)."""
    strands = _read(results_dir, image_id, "strands")
    strands = strands[strands["nucleus_id"] == int(nucleus_id)].sort_values("strand_id")
    p = Path(results_dir) / f"{image_id}{S.SUFFIXES['strips_' + kind]}"
    arr, desc = read_strips(p)
    ri = desc["roles"].index(role) if role in desc["roles"] else 0
    n_u = arr.shape[1]
    max_len = int(strands["strip_len"].max()) if len(strands) else 1
    gal = np.full((len(strands), n_u, max_len), np.nan, np.float32)
    for k, (_, r) in enumerate(strands.iterrows()):
        o, n = int(r["strip_offset"]), int(r["strip_len"])
        block = arr[ri, :, o:o + n].astype(np.float32)
        block[block == 0] = np.nan
        gal[k, :, :n] = block
    return gal, strands, desc


def strip_layers(results_dir, image_id: str, nucleus_id: int, kind: str = "flat") -> list[tuple]:
    gal, strands, desc = strip_gallery(results_dir, image_id, nucleus_id, kind)
    du = float(desc["du_um"] if kind == "flat" else desc["dv_um"])
    ds = float(desc["ds_um"])
    # rows stacked in y with a one-strip spacing; each strand's strip: (u, s)
    n, n_u, m = gal.shape
    canvas = np.full((max(n, 1) * (n_u + 2), m), np.nan, np.float32)
    for k in range(n):
        canvas[k * (n_u + 2):k * (n_u + 2) + n_u] = gal[k]
    layers = [(canvas, {"name": f"{image_id} nucleus {nucleus_id} strips ({kind})", "scale": (du, ds),
                        "translate": (0.0, float(desc["s_offset_um"])), "colormap": "magma"}, "image")]
    foci = _read(results_dir, image_id, "foci")
    if len(foci) and len(strands):
        f = foci[(foci["nucleus_id"] == int(nucleus_id)) & foci["assignment"].isin(["unique", "ambiguous"])]
        rows = {int(sid): k for k, sid in enumerate(strands["strand_id"])}
        pts = []
        for _, r in f.iterrows():
            k = rows.get(int(r["strand_id"])) if pd.notna(r["strand_id"]) else None
            if k is None:
                continue
            u = float(r["u_um"]) if kind == "flat" else float(r["v_um"])
            pts.append((k * (n_u + 2) * du + (n_u // 2) * du + u, float(r["s_um"])))
        if pts:
            layers.append((np.asarray(pts), {"name": f"{image_id} nucleus {nucleus_id} foci on strips", "size": 0.12,
                                             "face_color": "#00ffff"}, "points"))
    return layers


def image_layers(crop: np.ndarray, origin_vox, spacing, name: str = "SYP") -> list[tuple]:
    sp = tuple(float(s) for s in spacing)
    tr = tuple(float(o) * s for o, s in zip(origin_vox, sp))
    return [(crop, {"name": name, "scale": sp, "translate": tr, "colormap": "gray", "blending": "additive"}, "image")]
