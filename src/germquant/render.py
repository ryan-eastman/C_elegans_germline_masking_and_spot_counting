"""QC montages (and the entry point for napari renders later).

A per-image montage = max-intensity projection of each channel + a DNA/segmentation
overlay with nucleus boundaries + a scale bar. These double as supplementary figures and
let you eyeball every unattended result. Heavy/3D renders go through napari (viz extra).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from .fsutil import long_path


def _mip(img: np.ndarray) -> np.ndarray:
    return img.max(axis=0) if img.ndim == 3 else img


def _norm(a: np.ndarray) -> np.ndarray:
    a = a.astype(np.float32)
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)


def make_montage(
    stack,
    labels: np.ndarray | None,
    role_to_idx: dict[str, int | None],
    out_path: str | Path,
    *,
    foci_df=None,
    excluded_ids=None,
    sc_mask=None,
    granule_mask=None,
    scalebar_um: float = 10.0,
    title: str = "",
) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    panels = [(r, i) for r, i in role_to_idx.items() if i is not None]
    show_coloc = sc_mask is not None and granule_mask is not None
    n = len(panels) + (1 if labels is not None else 0) + (1 if show_coloc else 0)
    n = max(n, 1)

    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4))
    if n == 1:
        axes = [axes]

    dy = stack.spacing[1]
    bar_px = scalebar_um / dy

    k = 0
    for role, idx in panels:
        ax = axes[k]
        k += 1
        mip = _norm(_mip(stack.data[idx]))
        ax.imshow(mip, cmap="gray")
        ax.set_title(f"{role} (ch{idx})", fontsize=9)
        _scalebar(ax, mip.shape, bar_px, scalebar_um)
        ax.axis("off")

    if labels is not None:
        ax = axes[k]
        dna_idx = role_to_idx.get("dna")
        base = _norm(_mip(stack.data[dna_idx])) if dna_idx is not None else np.zeros(stack.shape_zyx[1:])
        ax.imshow(base, cmap="gray")
        lab_mip = labels.max(axis=0)
        if excluded_ids:
            # green = germline (kept), red = off-germline (dropped). inner boundaries carry the label.
            bnd = find_boundaries(lab_mip, mode="inner")
            is_excl = np.isin(lab_mip, list(excluded_ids))
            overlay = np.zeros((*bnd.shape, 4))
            overlay[bnd & ~is_excl] = (0, 1, 0, 1)
            overlay[bnd & is_excl] = (1, 0, 0, 1)
        else:
            bnd = find_boundaries(lab_mip, mode="outer")
            overlay = np.zeros((*bnd.shape, 4))
            overlay[bnd] = (1, 1, 0, 1)  # yellow nucleus outlines
        ax.imshow(overlay)
        if foci_df is not None and len(foci_df):
            ax.scatter(
                foci_df["x_um"] / stack.spacing[2], foci_df["y_um"] / stack.spacing[1],
                s=6, facecolors="none", edgecolors="magenta", linewidths=0.5,
            )
        n_excl = len(excluded_ids) if excluded_ids else 0
        seg_t = (f"germline={int(labels.max()) - n_excl} (red=off-gonad {n_excl})" if excluded_ids
                 else f"nuclei={int(labels.max())}")
        ax.set_title(seg_t + (f"  foci={len(foci_df)}" if foci_df is not None else ""), fontsize=9)
        _scalebar(ax, base.shape, bar_px, scalebar_um)
        ax.axis("off")

    if show_coloc:
        ax = axes[k]
        k += 1
        syp_idx = role_to_idx.get("central_element")
        base = (_norm(_mip(stack.data[syp_idx])) if syp_idx is not None
                else np.zeros(stack.shape_zyx[1:]))
        ax.imshow(base, cmap="gray")
        sc_mip = np.asarray(sc_mask).max(axis=0) > 0
        gr_mip = np.asarray(granule_mask).max(axis=0) > 0
        overlay = np.zeros((*sc_mip.shape, 4))
        overlay[sc_mip] = (1, 0, 0, 0.6)            # SYP aggregate = red
        overlay[gr_mip] = (0, 1, 0, 0.6)            # PGL-1 granules = green
        overlay[sc_mip & gr_mip] = (1, 1, 0, 0.95)  # overlap = yellow
        ax.imshow(overlay)
        ax.set_title(f"coloc: SYP-agg∩PGL-1 ({int((sc_mip & gr_mip).sum())} px)", fontsize=9)
        _scalebar(ax, base.shape, bar_px, scalebar_um)
        ax.axis("off")

    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    out_path = Path(out_path)
    os.makedirs(long_path(out_path.parent), exist_ok=True)
    fig.savefig(long_path(out_path), dpi=140, bbox_inches="tight", format="png")
    plt.close(fig)
    return out_path


def _scalebar(ax, shape, bar_px, um):
    h, w = shape
    y = h - max(8, h * 0.04)
    x0 = w * 0.05
    ax.plot([x0, x0 + bar_px], [y, y], "-", color="white", lw=3)
    ax.text(x0, y - h * 0.02, f"{um:g} µm", color="white", fontsize=8)
