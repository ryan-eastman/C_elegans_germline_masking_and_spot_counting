"""Grant pipeline figure as an editable SVG (Illustrator): artboard 180 x 48 mm (full text width, about
1/5 of a US-letter page), five steps left to right, real crops embedded as PNG, every element grouped
and labelled. Spare parts (alternative crops, single channels, whole-gonad views, arrows, legends,
formula, outputs box, alternative titles) sit OUTSIDE the artboard so they can be swapped in.
Run fig_grant_assets.py first. Output: coloc_analysis/grant/pipeline_figure.svg (+ copy on Desktop)."""
import base64
import os
import shutil

from PIL import Image

A = r"C:/Users/ryane/coloc_analysis/grant/assets"
OUT = r"C:/Users/ryane/coloc_analysis/grant/pipeline_figure.svg"
FONT = "Arial, Helvetica, sans-serif"
GREY = "#5A6270"
LIGHT = "#C9CED6"
ARROW = "#9AA3AD"
RED = "#E64A3C"
ZONE = {"early": "#56B4E9", "mid": "#009E73", "late": "#CC79A7"}
AW, AH = 180.0, 48.0          # artboard, mm
BOXW, GAP = 31.0, 5.75
IMG = 20.0                    # image edge, mm
X0 = [2 + i * (BOXW + GAP) for i in range(5)]
Y_TITLE, Y_IMG, Y_CAP = 4.6, 6.2, 29.6
T_TITLE, T_CAP, T_SMALL = 2.45, 1.95, 1.7   # font sizes in mm (2.45 mm ~ 7 pt, 1.95 mm ~ 5.5 pt)


def img_uri(name):
    with open(os.path.join(A, name), "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()


def img_size(name):
    with Image.open(os.path.join(A, name)) as im:
        return im.size


def image(name, x, y, h, gid, max_w=None):
    """embed PNG scaled to height h (mm), keeping aspect; centred if max_w given."""
    w0, h0 = img_size(name)
    w = h * w0 / h0
    if max_w is not None and w > max_w:
        w, h = max_w, max_w * h0 / w0
    xo = x + ((max_w - w) / 2 if max_w is not None else 0)
    return (f'<image id="{gid}" x="{xo:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'preserveAspectRatio="none" xlink:href="{img_uri(name)}"/>'), w, h


def text(x, y, s, size=T_CAP, weight="normal", anchor="start", fill="#111", style=""):
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" fill="{fill}" {style}>{s}</text>')


def lines(x, y, rows, size=T_CAP, dy=None, anchor="start", fill="#222"):
    dy = dy or size * 1.28
    return "\n".join(text(x, y + i * dy, r, size, anchor=anchor, fill=fill) for i, r in enumerate(rows))


def chevron(x, y, w=5.0, h=5.0, fill=ARROW, gid=""):
    return (f'<polygon id="{gid}" fill="{fill}" points="{x:.2f},{y - h / 2:.2f} {x + w * 0.62:.2f},{y - h / 2:.2f} '
            f'{x + w:.2f},{y:.2f} {x + w * 0.62:.2f},{y + h / 2:.2f} {x:.2f},{y + h / 2:.2f} {x + w * 0.3:.2f},{y:.2f}"/>')


def scalebar(x, y, um, px_per_mm_um, label=True):
    """um scale bar for a 24-um image drawn IMG mm wide: mm per um = IMG/24."""
    L = um * px_per_mm_um
    s = f'<line x1="{x:.2f}" y1="{y:.2f}" x2="{x + L:.2f}" y2="{y:.2f}" stroke="#fff" stroke-width="0.35"/>'
    if label:
        s += text(x + L / 2, y - 0.6, f"{um} µm", T_SMALL, anchor="middle", fill="#fff")
    return s


parts = []
P = parts.append
P('<?xml version="1.0" encoding="UTF-8"?>')
P(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
  f'width="275mm" height="165mm" viewBox="0 0 275 165" font-family="{FONT}">')
P('<defs><marker id="arrowhead" markerWidth="4" markerHeight="4" refX="3.2" refY="2" orient="auto">'
  f'<path d="M0,0 L4,2 L0,4 z" fill="{GREY}"/></marker></defs>')
P('<rect width="275" height="165" fill="#fff"/>')

# ---- artboard guide (delete in Illustrator) ----------------------------------------------------
P('<g id="artboard_guide">')
P(f'<rect x="0" y="0" width="{AW}" height="{AH}" fill="none" stroke="{LIGHT}" stroke-width="0.25" stroke-dasharray="1.5,1"/>')
P(text(0, AH + 3.2, f"artboard guide: {AW:.0f} x {AH:.0f} mm = full text width x 1/5 of a US-letter page (delete this guide)", T_SMALL, fill=GREY))
P('</g>')

# ---- steps ----------------------------------------------------------------------------------------
mm_per_um = IMG / 24.0
steps = [
    ("step1_input", "1  Confocal z-stack", "merge_24um.png",
     ["DAPI, PGL-1::GFP, SYP-3::mCherry,", "LMN-1; spinning disk, 0.2 µm z-steps", "ccw77 male and hermaphrodite,", "no HS vs heat shock"]),
    ("step2_nuclei", "2  Nuclear masks", "lamin_masks_24um.png",
     ["Cellpose on DAPI, then watershed", "to the LMN-1 envelope (red)", "nuclear interior excluded", "from all measurements"]),
    ("step3_staging", "3  Pachytene staging", "gonad_zones.png",
     ["hand-traced gonad axis", "start/end by the lab's row rule", "thirds: early, mid, late"]),
    ("step4_granules", "4  P granules", "granules_24um.png",
     ["Imaris-calibrated segmentation", "of PGL-1 (yellow) within a", "2.5 µm perinuclear shell (cyan)"]),
]
for i, (gid, title, png, cap) in enumerate(steps):
    x = X0[i]
    P(f'<g id="{gid}">')
    P(text(x, Y_TITLE, title, T_TITLE, "bold"))
    if png == "gonad_zones.png":
        im, w, h = image(png, x, Y_IMG, IMG, gid + "_img", max_w=BOXW)
        P(im)
        # zone legend swatches under the image
        lx = x
        for k, (zn, col) in enumerate(ZONE.items()):
            P(f'<rect x="{lx + k * 9.6:.2f}" y="{Y_CAP + 3.6:.2f}" width="2" height="2" fill="{col}"/>')
            P(text(lx + k * 9.6 + 2.6, Y_CAP + 5.35, zn, T_SMALL))
        P(lines(x, Y_CAP, cap[:2]))
        P(text(x, Y_CAP + 8.4, cap[2]))
    else:
        im, w, h = image(png, x, Y_IMG, IMG, gid + "_img")
        P(im)
        P(scalebar(x + 1.2, Y_IMG + IMG - 1.4, 5, mm_per_um))
        P(lines(x, Y_CAP, cap))
    P('</g>')
    P(chevron(x + BOXW + 0.4, Y_IMG + IMG / 2, w=GAP - 0.8, gid=f"arrow_{i + 1}"))

# step 5: vector schematic of the partition coefficient
x = X0[4]
P('<g id="step5_partition">')
P(text(x, Y_TITLE, "5  Partition coefficient", T_TITLE, "bold"))
cx, cy = x + IMG / 2, Y_IMG + IMG / 2
P(f'<rect x="{x:.2f}" y="{Y_IMG:.2f}" width="{IMG}" height="{IMG}" fill="#0B0B0B"/>')
P(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="8.2" fill="none" stroke="#19D3E6" stroke-width="0.3" stroke-dasharray="0.8,0.5"/>')   # shell edge 2.5 um
for r in (6.2, 7.4):                                                                                                                  # distance shells
    P(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r}" fill="none" stroke="#19D3E6" stroke-width="0.15" stroke-dasharray="0.4,0.5" opacity="0.6"/>')
P(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="5.0" fill="#2A2A2A" stroke="{RED}" stroke-width="0.45"/>')                               # nucleus, envelope
P(text(cx, cy + 0.7, "nucleus", T_SMALL, anchor="middle", fill="#BBB"))
import math
for k, (ang, rr) in enumerate([(20, 6.6), (75, 7.9), (140, 6.9), (200, 8.0), (250, 6.5), (300, 7.6)]):
    gx, gy = cx + rr * math.cos(math.radians(ang)), cy + rr * math.sin(math.radians(ang))
    P(f'<circle cx="{gx:.2f}" cy="{gy:.2f}" r="0.55" fill="#FFE11A"/>')
for ang, rr in [(45, 6.6), (110, 7.9), (170, 6.9), (225, 8.0), (275, 6.5), (330, 7.6)]:
    gx, gy = cx + rr * math.cos(math.radians(ang)), cy + rr * math.sin(math.radians(ang))
    P(f'<circle cx="{gx:.2f}" cy="{gy:.2f}" r="0.55" fill="none" stroke="#CFCFCF" stroke-width="0.2"/>')
P(text(x + 0.6, Y_IMG + IMG - 0.9, "2.5 µm shell, 0.25 µm bins", T_SMALL, fill="#19D3E6"))
# formula and controls
P(text(x, Y_CAP, "PC =", T_CAP, "bold"))
P(text(x + 5.2, Y_CAP - 1.05, "granule SYP-3 - bkg", 1.6))
P(f'<line x1="{x + 5.2:.2f}" y1="{Y_CAP - 0.35:.2f}" x2="{x + 27.5:.2f}" y2="{Y_CAP - 0.35:.2f}" stroke="#111" stroke-width="0.2"/>')
P(text(x + 5.2, Y_CAP + 1.35, "matched cytoplasm SYP-3 - bkg", 1.6))
P(lines(x, Y_CAP + 4.2, ["yellow: granule; open: cytoplasm at", "the same distance from the envelope", "controls: rotation null, z-shift floor", "output: PC per gonad and per zone"]))
P('</g>')

# ---- spares (outside the artboard) ------------------------------------------------------------------
P('<g id="spares">')
P(text(190, 4, "SPARE PARTS (outside the artboard): swap in / delete", T_TITLE, "bold", fill=GREY))
sp = [("merge_40um.png", "merged, 40 µm field"), ("lamin_masks_40um.png", "LMN-1 + masks, 40 µm"),
      ("masks_dapi_vs_lamin_24um.png", "DAPI label vs envelope"), ("granules_on_merge_24um.png", "granules on merge"),
      ("dapi_24um.png", "DAPI"), ("pgl_24um.png", "PGL-1::GFP"), ("syp_24um.png", "SYP-3::mCherry"), ("lamin_raw_24um.png", "LMN-1, no overlay"),
      ("merge_24um.png", "merged (copy)"), ("lamin_masks_24um.png", "LMN-1 + masks (copy)"), ("granules_24um.png", "granules (copy)")]
for k, (png, lab) in enumerate(sp):
    col, row = k % 4, k // 4
    sx, sy = 190 + col * 21.5, 6 + row * 25
    im, w, h = image(png, sx, sy, 18, f"spare_{os.path.splitext(png)[0]}")
    P(im)
    P(text(sx, sy + 18 + 2.2, lab, T_SMALL, fill=GREY))
# whole-gonad spares, bottom row
for k, (png, lab) in enumerate([("gonad_merge.png", "whole gonad, merged max projection"),
                                ("gonad_lamin_masks.png", "whole gonad, LMN-1 + envelope masks"),
                                ("gonad_zones.png", "whole gonad, zones + trace (copy)")]):
    sx = 2 + k * 40
    im, w, h = image(png, sx, 62, 42, f"spare_{os.path.splitext(png)[0]}", max_w=36)
    P(im)
    P(text(sx, 62 + 42 + 2.4, lab, T_SMALL, fill=GREY))
# vector spares
P(text(130, 64, "arrows / boxes", T_CAP, "bold", fill=GREY))
P(chevron(130, 70, w=8, h=6, gid="spare_chevron_large"))
P(chevron(140, 70, w=5, h=5, gid="spare_chevron_small"))
P(f'<line x1="148" y1="70" x2="160" y2="70" stroke="{GREY}" stroke-width="0.4" marker-end="url(#arrowhead)" id="spare_arrow_line"/>')
P(f'<path d="M163,73 Q168,60 175,73" fill="none" stroke="{GREY}" stroke-width="0.4" marker-end="url(#arrowhead)" id="spare_arrow_curved"/>')
P(f'<rect x="130" y="76" width="30" height="14" rx="1.2" fill="#F3F5F8" stroke="{LIGHT}" stroke-width="0.3" id="spare_box"/>')
P(text(145, 84, "step box", T_CAP, anchor="middle", fill=GREY))
# zone legend spare
P('<g id="spare_zone_legend">')
for k, (zn, col) in enumerate(ZONE.items()):
    P(f'<rect x="{130 + k * 17:.2f}" y="94" width="2.2" height="2.2" fill="{col}"/>')
    P(text(130 + k * 17 + 2.8, 95.9, zn + " pachytene", T_SMALL))
P('</g>')
# channel legend spare
P('<g id="spare_channel_legend">')
for k, (nm, col) in enumerate([("DAPI", "#4C6FFF"), ("PGL-1::GFP", "#2ECC40"), ("SYP-3::mCherry", "#FF4136"), ("LMN-1", "#DDDDDD")]):
    P(f'<rect x="{130 + k * 19:.2f}" y="100" width="2.2" height="2.2" fill="{col}" stroke="#888" stroke-width="0.1"/>')
    P(text(130 + k * 19 + 2.8, 101.9, nm, T_SMALL))
P('</g>')
# outputs box spare
P('<g id="spare_outputs_box">')
P(f'<rect x="130" y="106" width="60" height="22" rx="1.2" fill="#fff" stroke="{GREY}" stroke-width="0.3"/>')
P(text(132, 110, "Outputs", T_CAP, "bold"))
P(lines(132, 113.5, ["partition coefficient per gonad and per pachytene third", "granule-specific enrichment (PC / z-shift floor)",
                     "fraction of granules with SYP-3 above local cytoplasm", "male vs hermaphrodite, no HS vs heat shock"], T_SMALL))
P('</g>')
# mini dot plot icon spare
P('<g id="spare_dotplot_icon">')
P(f'<rect x="195" y="106" width="22" height="22" fill="#fff" stroke="{LIGHT}" stroke-width="0.25"/>')
P('<line x1="199" y1="125" x2="215" y2="125" stroke="#333" stroke-width="0.25"/><line x1="199" y1="109" x2="199" y2="125" stroke="#333" stroke-width="0.25"/>')
for gx, ys, col in [(203.5, [121, 122.5, 120.2, 121.8], "#0072B2"), (210.5, [114, 116.5, 112.5, 115.2, 117.5], "#D55E00")]:
    for yv in ys:
        P(f'<circle cx="{gx + (yv * 7) % 3 - 1.5:.2f}" cy="{yv:.2f}" r="0.55" fill="none" stroke="{col}" stroke-width="0.25"/>')
    P(f'<line x1="{gx - 2.2:.2f}" y1="{sum(ys) / len(ys):.2f}" x2="{gx + 2.2:.2f}" y2="{sum(ys) / len(ys):.2f}" stroke="#111" stroke-width="0.35"/>')
P(text(203.5, 127.5, "no HS", T_SMALL, anchor="middle")); P(text(210.5, 127.5, "HS", T_SMALL, anchor="middle"))
P('</g>')
# alternative titles / short captions spare
P('<g id="spare_alt_text">')
P(text(225, 108, "alternative titles / captions", T_CAP, "bold", fill=GREY))
P(lines(225, 111.5, ["Image acquisition", "Lamin-envelope nuclear masks", "Staging by hand-traced axis", "Perinuclear P-granule segmentation",
                     "Distance-matched partition coefficient", "Rotation null and z-shift controls", "Per-gonad and per-zone quantification",
                     "n = 13 gonads (9 male, 4 hermaphrodite)"], T_SMALL, dy=2.6))
P('</g>')
# formula spare (large)
P('<g id="spare_formula">')
P(text(2, 112, "PC =", 3.0, "bold"))
P(text(9.5, 110.4, "SYP-3 in granule - background", T_CAP))
P('<line x1="9.5" y1="111.3" x2="60" y2="111.3" stroke="#111" stroke-width="0.25"/>')
P(text(9.5, 114.2, "SYP-3 in cytoplasm at the same distance from the envelope - background", T_CAP))
P(text(2, 119, "granule-specific enrichment = PC / PC of the same masks shifted 1.2 µm in z (no granule under the mask)", T_SMALL))
P('</g>')
P('</g>')   # spares
P('</svg>')

os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, "w", encoding="utf-8").write("\n".join(parts))
print("wrote", OUT, f"{os.path.getsize(OUT) // 1024} KB")
for d in [r"C:/Users/ryane/OneDrive/Desktop/germquant_figures"]:      # same folder pubstyle.save uses
    if os.path.isdir(d):
        shutil.copy(OUT, os.path.join(d, "pipeline_figure_grant.svg"))
        print("copied to", d)
