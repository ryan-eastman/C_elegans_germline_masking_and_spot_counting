"""Foci on the straightened strands: each spot of a spot table (RAD-51 ``spots``, COSA-1 or another
``spots_<name>``) is mapped to the closest point of the strands of its own nucleus.

A spot is ``unique`` when |u| <= ``foci.max_um[1]`` (lateral) and |v| <= ``foci.max_um[0]`` (axial) and it lies
within 0.3 um of the strand's ends along the tangent; ``ambiguous`` when a second strand's PSF-scaled distance is
within ``foci.ambiguity_ratio`` of the first's; otherwise ``off_strand`` (``no_strand`` when the nucleus has no
strand). Both PSF-scaled distances count the lateral (u), axial (v) and along-tangent offsets from the strand
proper (s in [0, L]), so a strand that ends short of a focus and only points at it is not a second candidate.
s = 0 is the end with the smaller rounded (z, y, x), so folded positions and inter-focus distances do not depend
on the orientation.
"""
from __future__ import annotations

import numpy as np

from .frames import beyond_end as _beyond
from .frames import to_strip

FOCI_FIELDS = ["spot_table", "spot_id", "nucleus_id", "marker", "assignment", "strand_key", "s_um", "s_frac",
               "s_folded_frac", "u_um", "v_um", "dist_um", "second_strand_key", "in_gap", "rank_on_strand",
               "next_focus_um", "next_focus_frac", "z_um", "y_um", "x_um"]


def assign_foci(spots: dict, frames: list, lines: list, gaps: list, psf, p, nid: int) -> list[dict]:
    max_ax, max_lat = (float(v) for v in p["foci.max_um"])
    ratio = float(p["foci.ambiguity_ratio"])
    sz, sxy = float(psf[0]), float(psf[1])
    rows: list[dict] = []
    for table, df in sorted(spots.items()):
        if df is None or not len(df):
            continue
        pts = df[["z_um", "y_um", "x_um"]].to_numpy(float)
        res = to_strip(frames, pts) if frames else None
        for i in range(len(df)):
            r = {"spot_table": table, "spot_id": df["spot_id"].iloc[i] if "spot_id" in df else i, "nucleus_id": nid,
                 "marker": df["marker"].iloc[i] if "marker" in df else "", "assignment": "no_strand",
                 "strand_key": -1, "s_um": np.nan, "s_frac": np.nan, "s_folded_frac": np.nan, "u_um": np.nan,
                 "v_um": np.nan, "dist_um": np.nan, "second_strand_key": -1, "in_gap": False,
                 "rank_on_strand": np.nan, "next_focus_um": np.nan, "next_focus_frac": np.nan,
                 "z_um": pts[i, 0], "y_um": pts[i, 1], "x_um": pts[i, 2]}
            if res is not None and res["strand"][i] >= 0:
                j = int(res["strand"][i])
                L = lines[j].L
                s, u, v, w = res["s"][i], res["u"][i], res["v"][i], res["w"][i]
                nz = float(abs(frames[j].frame_at(s)[3][0][0]))
                sax = float(np.sqrt(sxy ** 2 * (1 - nz ** 2) + sz ** 2 * nz ** 2))
                d1 = float(np.sqrt((u / sxy) ** 2 + (v / sax) ** 2 + (_beyond(s, w, L) / sxy) ** 2))
                r.update({"strand_key": j, "s_um": s, "s_frac": float(np.clip(s / L, 0, 1)),
                          "s_folded_frac": float(min(max(s, 0.0), max(L - s, 0.0)) / L), "u_um": u, "v_um": v,
                          "dist_um": res["dist"][i], "second_strand_key": int(res["second"][i]),
                          "in_gap": bool(any(g["s_start_um"] <= s <= g["s_end_um"] for g in gaps[j]))})
                on = abs(u) <= max_lat and abs(v) <= max_ax and (-0.3 <= s <= L + 0.3) and abs(w) <= 0.3
                if not on:
                    r["assignment"] = "off_strand"
                else:
                    amb = False
                    k = int(res["second"][i])
                    if k >= 0:
                        # the second strand's distance on the same PSF scale (its own frame is near-identical
                        # at this scale), including the distance beyond its end along its tangent: a strand
                        # that only points at the focus is not a second candidate
                        along2 = _beyond(res["second_s"][i], res["second_w"][i], lines[k].L)
                        d2 = float(np.sqrt((res["second_u"][i] / sxy) ** 2 + (res["second_v"][i] / sax) ** 2
                                           + (along2 / sxy) ** 2))
                        amb = np.isfinite(d2) and d2 <= ratio * max(d1, 1e-9)
                    r["assignment"] = "ambiguous" if amb else "unique"
            rows.append(r)
    # ranks and spacing along each strand (assigned foci, per table)
    by = {}
    for r in rows:
        if r["assignment"] in ("unique", "ambiguous"):
            by.setdefault((r["spot_table"], r["strand_key"]), []).append(r)
    for (_t, j), group in by.items():
        group.sort(key=lambda r: (r["s_um"], str(r["spot_id"])))
        L = lines[j].L
        for k, r in enumerate(group):
            r["rank_on_strand"] = k + 1
            if k + 1 < len(group):
                nxt = float(group[k + 1]["s_um"] - r["s_um"])
                r["next_focus_um"] = nxt
                r["next_focus_frac"] = nxt / L if L > 0 else np.nan
    return rows
