"""Pop-up tool: draw the pachytene region on each gonad with a multi-point line.

Ryan's design (2026-08-23): automated pachytene-start detection plateaued at ~5 rows error, so the
human draws the region instead - click a polyline along the gonad from PACHYTENE START to PACHYTENE END
(as many points as the curvature needs). Everything downstream is automated:
  * every germline nucleus is projected onto the drawn line -> arc-length position + perpendicular dist
  * zones = equal-length thirds of the drawn line (early / mid / late), the lab definition
  * nuclei before the first point = "pre", beyond the last = "post", farther than OFF_AXIS_UM from the
    line = "off_axis" (this excludes fused second arms / debris BY CONSTRUCTION - the failure mode no
    automated axis could solve)

USAGE
  .venv\Scripts\python.exe trace_pachytene.py            # loop over every staged gonad not yet traced
  .venv\Scripts\python.exe trace_pachytene.py <image_id> # one specific gonad
  .venv\Scripts\python.exe trace_pachytene.py --redo     # loop over ALL, including already-traced

CONTROLS (shown in the window title too)
  left-click   add a point (start clicking at pachytene START, distal side)
  z            undo last point
  r            reset this gonad
  enter / n    accept trace -> saves + next gonad
  s            skip this gonad (marked unusable)
  b            back to previous gonad
  d            toggle SYP-3 (red) overlay
  c            toggle nucleus centroids
  q            quit (progress is saved after every gonad)
  (use the toolbar magnifier to zoom; clicks are ignored while zoom/pan is active)

OUTPUT (all in coloc_analysis/staging/)
  pachytene_traces.json          polyline per gonad, crop-um coordinates, resumable
  zones/<iid>_zones.csv          per-nucleus: s_um along the line, r_um off it, zone
Coordinate frame: the same crop-um frame as staging/<iid>_nuclei.csv (cx, cy), i.e. the cached crops."""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")

STAGING = r"C:/Users/ryane/coloc_analysis/staging"
TRACES = os.path.join(STAGING, "pachytene_traces.json")
ZONE_DIR = os.path.join(STAGING, "zones")
OFF_AXIS_UM = 20.0        # nuclei farther than this from the drawn line are excluded
PXY = 0.1083              # um / px in xy
DS = 2                    # display downsample
os.makedirs(ZONE_DIR, exist_ok=True)


# ----------------------------- geometry (headless-testable) -----------------------------
def project_to_polyline(pts, poly):
    """pts (n,2) and poly (k,2) in um -> (s, r): arc-length along poly (unclamped at the ends,
    so s<0 means before the start and s>L means beyond the end) and perpendicular distance."""
    pts = np.asarray(pts, float)
    poly = np.asarray(poly, float)
    seg = np.diff(poly, axis=0)
    seglen = np.linalg.norm(seg, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seglen)])
    best_r = np.full(len(pts), np.inf)
    best_s = np.zeros(len(pts))
    for i in range(len(seg)):
        d = seg[i]
        L2 = max(seglen[i] ** 2, 1e-12)
        t = ((pts - poly[i]) @ d) / L2
        tc = np.clip(t, 0.0, 1.0)
        proj = poly[i] + tc[:, None] * d
        r = np.linalg.norm(pts - proj, axis=1)
        upd = r < best_r
        s = cum[i] + tc * seglen[i]
        # unclamped overshoot at the free ends so pre/post can be distinguished
        if i == 0:
            s = np.where(t < 0, cum[i] + t * seglen[i], s)
        if i == len(seg) - 1:
            s = np.where(t > 1, cum[i] + t * seglen[i], s)
        best_s = np.where(upd, s, best_s)
        best_r = np.where(upd, r, best_r)
    return best_s, best_r, float(cum[-1])


def assign_zones(df, poly, off_axis_um=OFF_AXIS_UM, adaptive=True):
    """df needs cx, cy (crop um). Returns df + s_um, r_um, zone, off_axis_cut.
    adaptive: the off-axis cutoff scales with the gonad's own tube width, cut = min(off_axis_um,
    2.5 x median r of nuclei inside the traced span). A fixed 20 um readmitted second-arm nuclei in
    narrow male gonads (noHS_male_003: 45 nuclei at 12-20 um) while herm tubes genuinely reach ~18 um."""
    s, r, L = project_to_polyline(df[["cx", "cy"]].to_numpy(), poly)
    cut = float(off_axis_um)
    if adaptive:
        inside0 = (r <= off_axis_um) & (s >= 0) & (s <= L)
        if inside0.sum() >= 20:
            cut = float(min(off_axis_um, 2.5 * np.median(r[inside0])))
    zone = np.full(len(df), "off_axis", dtype=object)
    on = r <= cut
    f = s / max(L, 1e-9)
    zone[on & (s < 0)] = "pre"
    zone[on & (s > L)] = "post"
    inside = on & (s >= 0) & (s <= L)
    zone[inside & (f < 1 / 3)] = "early"
    zone[inside & (f >= 1 / 3) & (f < 2 / 3)] = "mid"
    zone[inside & (f >= 2 / 3)] = "late"
    out = df.copy()
    out["s_um"], out["r_um"], out["zone"] = s, r, zone
    out["off_axis_cut_um"] = round(cut, 2)
    return out, L


# ----------------------------- persistence -----------------------------
def load_traces():
    if os.path.exists(TRACES):
        return json.load(open(TRACES, encoding="utf-8"))
    return {}


def save_trace(iid, points, status):
    tr = load_traces()
    tr[iid] = {"points_um": [[round(float(x), 2), round(float(y), 2)] for x, y in points],
               "status": status, "n_points": len(points),
               "off_axis_um": OFF_AXIS_UM, "traced_at": time.strftime("%Y-%m-%d %H:%M")}
    with open(TRACES, "w", encoding="utf-8") as f:
        json.dump(tr, f, indent=1)
    if status == "skipped":
        # retire any zones from an earlier trace of this gonad so no worker can consume them
        zf = os.path.join(ZONE_DIR, f"{iid}_zones.csv")
        if os.path.exists(zf):
            os.replace(zf, zf + ".retired")
    if status == "traced" and len(points) >= 2:
        df = pd.read_csv(os.path.join(STAGING, f"{iid}_nuclei.csv"))
        zoned, L = assign_zones(df, np.asarray(points))
        zoned.to_csv(os.path.join(ZONE_DIR, f"{iid}_zones.csv"), index=False)
        n = zoned.zone.value_counts().to_dict()
        print(f"  saved {iid[-24:]}: L={L:.0f} um, zones={n}")


# ----------------------------- the pop-up -----------------------------
def run_gui(iids):
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

    from skimage import exposure

    state = {"idx": 0, "pts": [], "show_syp": True, "show_cent": False,
             "zmode": "mid", "zi": 0, "clahe": True, "gamma": 1.0}
    fig, ax = plt.subplots(figsize=(15, 10))
    try:
        fig.canvas.manager.set_window_title("Pachytene tracer")
    except Exception:
        pass

    cache = {"iid": None, "stack": None, "views": {}}

    def mips(iid):
        """(dapi, syp) MIP for the current z-mode, CLAHE-enhanced DAPI. z-slabs stop overlapping
        cell layers from hiding each other; CLAHE stops the bright distal mass from blowing out."""
        if cache["iid"] != iid:
            z = np.load(os.path.join(STAGING, "cache", f"{iid}.npz"))
            cache.update(iid=iid, stack=(z["dapi"][:, ::DS, ::DS], z["syp"][:, ::DS, ::DS]), views={})
            state["zi"] = z["dapi"].shape[0] // 2
        D, S = cache["stack"]
        nz = D.shape[0]
        plane = state["zmode"] == "plane"
        key = ("plane", state["zi"], state["clahe"]) if plane else (state["zmode"], state["clahe"])
        if not plane and key in cache["views"]:
            return cache["views"][key]
        if plane:
            zi = int(np.clip(state["zi"], 0, nz - 1))
            d = D[zi].astype(np.float32)
            s = S[zi].astype(np.float32)
        else:
            zsl = {"all": slice(None), "low": slice(0, max(1, nz // 3)),
                   "mid": slice(nz // 4, max(nz // 4 + 1, 3 * nz // 4)),
                   "high": slice(2 * nz // 3, nz)}[state["zmode"]]
            d = D[zsl].max(0).astype(np.float32)
            s = S[zsl].max(0).astype(np.float32)

        def st(a, lo, hi):
            p0, p1 = np.percentile(a, [lo, hi])
            return np.clip((a - p0) / (p1 - p0 + 1e-9), 0, 1)
        dn = st(d, 1.0, 99.7)
        if state["clahe"]:
            dn = exposure.equalize_adapthist(dn, kernel_size=128, clip_limit=0.01).astype(np.float32)
        view = (dn, st(s, 40.0, 99.7))
        if not plane:                      # single planes are cheap to recompute; caching them bloats RAM
            cache["views"][key] = view
        return view

    def draw():
        iid = iids[state["idx"]]
        keep_view = state.get("cur_iid") == iid and ax.images
        if keep_view:
            xl, yl = ax.get_xlim(), ax.get_ylim()
        state["cur_iid"] = iid
        dn, sn = mips(iid)
        dg = np.power(dn, state["gamma"]) if state["gamma"] != 1.0 else dn
        rgb = np.stack([dg, dg, dg], -1)
        if state["show_syp"]:
            rgb[..., 0] = np.maximum(rgb[..., 0], sn)
        ax.clear()
        h, w = dn.shape
        ax.imshow(rgb, extent=[0, w * PXY * DS, h * PXY * DS, 0], interpolation="bilinear")
        if keep_view:
            ax.set_xlim(xl), ax.set_ylim(yl)
        if state["show_cent"]:
            nc = pd.read_csv(os.path.join(STAGING, f"{iid}_nuclei.csv"))
            ax.scatter(nc.cx, nc.cy, s=4, c="#3af", alpha=0.5)
        P = state["pts"]
        if P:
            xs, ys = zip(*P)
            ax.plot(xs, ys, "-o", color="lime", lw=2, ms=5)
            ax.annotate("START", P[0], color="lime", fontsize=11, fontweight="bold",
                        xytext=(8, -8), textcoords="offset points")
            if len(P) > 1:
                ax.annotate("END", P[-1], color="magenta", fontsize=11, fontweight="bold",
                            xytext=(8, -8), textcoords="offset points")
        done = sum(1 for g in load_traces().values() if g["status"] in ("traced", "skipped"))
        nz = cache["stack"][0].shape[0] if cache["stack"] is not None else 0
        vdesc = f"plane {state['zi'] + 1}/{nz}" if state["zmode"] == "plane" else state["zmode"]
        ax.set_title(f"[{state['idx'] + 1}/{len(iids)}  traced:{done}]  {iid}    "
                     f"view: z={vdesc}  clahe={'on' if state['clahe'] else 'off'}  gamma={state['gamma']:.2f}\n"
                     "click: add point (start at pachytene START, end at pachytene END)   z:undo  r:reset  enter:accept  s:skip  b:back  q:quit\n"
                     "view:  SCROLL WHEEL = step through single z-planes   m or 4 = max projection   1/2/3 = lower/mid/upper z-slab   "
                     "v = CLAHE   up/down = brightness   d = SYP   c = centroids",
                     fontsize=9)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        fig.canvas.draw_idle()

    def toolbar_active():
        tb = getattr(fig.canvas.manager, "toolbar", None)
        return bool(getattr(tb, "mode", ""))

    def on_click(ev):
        if ev.inaxes != ax or toolbar_active() or ev.button != 1:
            return
        state["pts"].append((float(ev.xdata), float(ev.ydata)))
        draw()

    def nxt():
        state["pts"] = []
        if state["idx"] + 1 < len(iids):
            state["idx"] += 1
            draw()
        else:
            print("all gonads done")
            plt.close(fig)

    def on_key(ev):
        k = ev.key
        iid = iids[state["idx"]]
        if k == "z" and state["pts"]:
            state["pts"].pop(); draw()
        elif k == "r":
            state["pts"] = []; draw()
        elif k in ("enter", "n"):
            if len(state["pts"]) < 2:
                print("  need at least 2 points (start + end)"); return
            save_trace(iid, state["pts"], "traced"); nxt()
        elif k == "s":
            save_trace(iid, state["pts"], "skipped"); print(f"  skipped {iid[-24:]}"); nxt()
        elif k == "b" and state["idx"] > 0:
            state["idx"] -= 1; state["pts"] = []; draw()
        elif k == "d":
            state["show_syp"] = not state["show_syp"]; draw()
        elif k == "c":
            state["show_cent"] = not state["show_cent"]; draw()
        elif k in ("1", "2", "3", "4", "m"):
            state["zmode"] = {"1": "low", "2": "mid", "3": "high", "4": "all", "m": "all"}[k]; draw()
        elif k == "v":
            state["clahe"] = not state["clahe"]; draw()
        elif k == "up":
            state["gamma"] = max(0.3, state["gamma"] * 0.8); draw()
        elif k == "down":
            state["gamma"] = min(3.0, state["gamma"] * 1.25); draw()
        elif k == "q":
            plt.close(fig)

    def on_scroll(ev):
        """mouse wheel steps through single z-planes (enters plane mode from any view)."""
        nz = cache["stack"][0].shape[0] if cache["stack"] is not None else 1
        if state["zmode"] != "plane":
            state["zmode"] = "plane"
        state["zi"] = int(np.clip(state["zi"] + (1 if ev.step > 0 else -1), 0, nz - 1))
        draw()

    fig.canvas.mpl_connect("button_press_event", on_click)
    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("scroll_event", on_scroll)
    draw()
    plt.show()


def main():
    args = [a for a in sys.argv[1:]]
    redo = "--redo" in args
    args = [a for a in args if not a.startswith("--")]
    have = sorted(f[:-11] for f in os.listdir(STAGING) if f.endswith("_nuclei.csv")
                  and os.path.exists(os.path.join(STAGING, "cache", f[:-11] + ".npz")))
    if args:
        iids = [a for a in args if a in have]
        missing = [a for a in args if a not in have]
        for m in missing:
            print(f"no staged data for {m} (run crescent_axis.py on it first)")
    else:
        done = {k for k, v in load_traces().items() if v["status"] in ("traced", "skipped")}
        iids = have if redo else [i for i in have if i not in done]
    if not iids:
        print("nothing to trace (all done - use --redo to retrace)")
        return
    print(f"{len(iids)} gonads to trace")
    run_gui(iids)


if __name__ == "__main__":
    main()
