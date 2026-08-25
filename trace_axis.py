"""Manual germline-axis tracer (pure DISTAL -> PROXIMAL). For each gonad, YOU click a polyline from the
DISTAL tip to the PROXIMAL end along the tube; EVERY germline nucleus is projected onto that line
(arc length -> axis 0..1). Nothing is dropped. Pachytene sub-staging is done LATER from the SC/SYP
signal along this axis, so you only need to draw the geometry. Raw click-lines are saved to
axis_traces.json so a re-projection never needs re-tracing.

RUN (from the repo root, with a GUI available):
    .venv\\Scripts\\python trace_axis.py                # every gonad
    .venv\\Scripts\\python trace_axis.py --redo         # only the ones you traced before (axis_traced.csv)
    .venv\\Scripts\\python trace_axis.py HS_male_009 HS_herm_14   # only ids ENDING with these

CONTROLS per gonad window:
    left-click     add a point (DISTAL tip -> ... -> PROXIMAL end, following the curve/folds)
    backspace/del  undo last point
    ENTER          accept -> project, go to next
    s              SKIP (keep the current axis for this gonad)
    q              quit (save what's done so far)
Dots are colored by the current axis (viridis, dark=distal) with D/P as a reference; gray dots are
detached pieces the auto-fit dropped (trace through them to include them). First click = distal (0),
last = proximal (1).
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SP = r"C:/Users/ryane/AppData/Local/Temp/claude/C--Users-ryane/0120e6bd-0a6e-4982-a789-021eec984f0f/scratchpad"
DIRS = [r"C:/Users/ryane/ccw77_fullres", r"C:/Users/ryane/ccw77_0622", r"C:/Users/ryane/n2_0625"]
AXCSV = os.path.join(SP, "axis_fixed.csv")

nfiles = {}
for root in DIRS:
    for nf in glob.glob(os.path.join(root, "**", "*__nuclei.csv"), recursive=True):
        nfiles[os.path.basename(nf).replace("__nuclei.csv", "")] = nf

ax_df = pd.read_csv(AXCSV)
cur = {iid: sub.set_index("nucleus_id")["axis_new"] for iid, sub in ax_df.groupby("image_id")}
args = [a for a in sys.argv[1:] if a != "--redo"]
targets = sorted(nfiles)
if "--redo" in sys.argv:                                    # only the gonads traced last time
    done = pd.read_csv(os.path.join(SP, "axis_traced.csv")).image_id.tolist()
    targets = [iid for iid in targets if iid in done]
elif args:
    targets = [iid for iid in targets if any(iid.endswith(s) for s in args)]


def project_to_polyline(x, y, poly):
    """Arc length (0..1) of each nucleus's nearest point on the distal->proximal polyline. No dropping:
    every germline nucleus gets an axis (pachytene bounding happens later from the SC signal)."""
    P = np.asarray(poly, float)
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    uu = np.linspace(0, seg[-1], 800)
    dense = np.c_[np.interp(uu, seg, P[:, 0]), np.interp(uu, seg, P[:, 1])]
    d2 = (x[:, None] - dense[None, :, 0]) ** 2 + (y[:, None] - dense[None, :, 1]) ** 2
    return uu[d2.argmin(1)] / uu[-1]


state = {}
polylines = {}
for iid in targets:
    g = pd.read_csv(nfiles[iid])
    gg = g[g.in_germline.astype(bool)]                 # ALL germline nuclei (incl. detached pieces)
    x, y = gg.centroid_x_um.to_numpy(), gg.centroid_y_um.to_numpy()
    if iid in cur:
        a = cur[iid].reindex(gg.nucleus_id).to_numpy() # current axis; NaN for pieces the auto-fit dropped
    else:
        a = np.full(len(gg), np.nan)                   # previously-excluded gonad: no current axis
    fin = np.isfinite(a)

    fig, ax = plt.subplots(figsize=(9, 8))
    if (~fin).any():
        ax.scatter(x[~fin], y[~fin], c="lightgray", s=16)   # dropped pieces you can now include
    ax.scatter(x[fin], y[fin], c=a[fin], cmap="viridis", s=16, vmin=0, vmax=1)
    for lab, m in [("D", fin & (a < 0.06)), ("P", fin & (a > 0.94))]:
        if m.any():
            ax.text(x[m].mean(), y[m].mean(), lab, color="w", fontsize=15, weight="bold", ha="center",
                    va="center", bbox=dict(boxstyle="circle", fc="k", ec="w"))
    ax.invert_yaxis(); ax.set_aspect("equal")
    ax.set_title(f"{iid}\nclick DISTAL -> PROXIMAL along tube | ENTER=accept  s=skip  q=quit  backspace=undo")
    pts = []
    line, = ax.plot([], [], "-o", color="red", lw=2, ms=6)
    decision = {"action": "skip"}

    def redraw():
        if pts:
            line.set_data([p[0] for p in pts], [p[1] for p in pts])
        else:
            line.set_data([], [])
        fig.canvas.draw_idle()

    def on_click(ev):
        if ev.inaxes is ax and ev.button == 1 and ev.xdata is not None:
            pts.append((ev.xdata, ev.ydata)); redraw()

    def on_key(ev):
        if ev.key == "enter" and len(pts) >= 2:
            decision["action"] = "accept"; plt.close(fig)
        elif ev.key in ("backspace", "delete") and pts:
            pts.pop(); redraw()
        elif ev.key == "s":
            decision["action"] = "skip"; plt.close(fig)
        elif ev.key == "q":
            decision["action"] = "quit"; plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_click)
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()

    if decision["action"] == "quit":
        break
    if decision["action"] == "accept":
        axis = project_to_polyline(x, y, pts)
        state[iid] = pd.DataFrame({"image_id": iid, "nucleus_id": gg.nucleus_id.to_numpy(),
                                   "axis_new": axis})
        polylines[iid] = [[float(px), float(py)] for px, py in pts]
        print(f"traced {iid}: {len(pts)} pts -> {len(gg)} nuclei (distal->proximal, none dropped)")

if state:
    keep = ax_df[~ax_df.image_id.isin(state)]
    out = pd.concat([keep] + list(state.values()), ignore_index=True)
    out.to_csv(AXCSV, index=False)
    # merge raw click-lines into axis_traces.json (so re-projection never needs re-tracing)
    jf = os.path.join(SP, "axis_traces.json")
    prev = json.load(open(jf)) if os.path.exists(jf) else {}
    prev.update(polylines)
    json.dump(prev, open(jf, "w"), indent=0)
    pd.DataFrame({"image_id": list(state)}).to_csv(os.path.join(SP, "axis_traced.csv"), index=False)
    print(f"\nupdated axis_fixed.csv + axis_traces.json: traced {len(state)} gonads")
else:
    print("no gonads traced (all skipped)")
