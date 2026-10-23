"""Outputs of sc-straighten: atomic writers, the per-image record, the sub-hash and resume, the strips TIFFs, the
overview PNG, and the ``germquant collect`` hook (``collect_tables``).

Every file is new and written atomically (a temporary file in the same folder, then ``os.replace``, through the
extended-length path); the record ``<id>__sc_straighten.json`` is written last, so its ``done: true`` means every
output of that run is on disk.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

from ...fsutil import long_path
from . import schema as S

log = logging.getLogger(__name__)

TOOLS = ("numpy", "scipy", "scikit-image", "networkx", "tifffile", "h5py", "nd2", "pandas")


def tool_versions() -> dict:
    from importlib import metadata

    out = {}
    for t in TOOLS:
        try:
            out[t] = metadata.version(t)
        except Exception:  # noqa: BLE001
            out[t] = "not-installed"
    return out


def _tmp(path: Path) -> Path:
    return path.with_name(f".{path.name}.tmp{os.getpid()}")


def _replace(tmp: Path, path: Path) -> None:
    os.replace(long_path(tmp), long_path(path))


def write_text(path: str | Path, text: str) -> None:
    path = Path(path)
    tmp = _tmp(path)
    with open(long_path(tmp), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    _replace(tmp, path)


def conform(df: pd.DataFrame, declared: list[str]) -> pd.DataFrame:
    df = df.copy()
    for c in declared:
        if c not in df.columns:
            df[c] = pd.NA
    return df[list(declared) + [c for c in df.columns if c not in declared]]


def write_table(df: pd.DataFrame, path: str | Path, declared: list[str], meta: dict | None = None) -> None:
    df = conform(df, declared)
    for k, v in (meta or {}).items():
        df[k] = v
    path = Path(path)
    tmp = _tmp(path)
    df.to_csv(long_path(tmp), index=False, float_format="%.6g", lineterminator="\n")
    _replace(tmp, path)


def write_strips(path: str | Path, arr: np.ndarray, desc: dict) -> None:
    """(roles, rows, columns) float counts -> uint16 (rounded; 0 = outside the image), zlib, JSON description."""
    import tifffile

    a = np.nan_to_num(np.asarray(arr, np.float64), nan=0.0)
    a = np.clip(np.rint(a), 0, 65535).astype(np.uint16)
    if a.shape[-1] == 0:
        a = np.zeros(a.shape[:-1] + (1,), np.uint16)
        desc = dict(desc, empty=True)
    path = Path(path)
    tmp = _tmp(path)
    tifffile.imwrite(long_path(tmp), a, compression="zlib", description=json.dumps(desc), metadata=None)
    _replace(tmp, path)


def read_strips(path: str | Path) -> tuple[np.ndarray, dict]:
    import tifffile

    with tifffile.TiffFile(long_path(path)) as tif:
        arr = tif.asarray()
        desc = json.loads(tif.pages[0].description)
    if arr.ndim == 2:
        arr = arr[None]
    return arr, desc


def read_record(path: str | Path) -> dict | None:
    try:
        with open(long_path(path), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, Path):
        return str(o)
    if o is pd.NA:
        return None
    return o


def write_record(path: str | Path, record: dict) -> None:
    write_text(path, json.dumps(_jsonable(record), indent=1, sort_keys=False))


def sub_hash(payload: dict) -> str:
    return hashlib.sha256(json.dumps(_jsonable(payload), sort_keys=True, separators=(",", ":"),
                                     default=str).encode()).hexdigest()[:12]


def paths_for(out_dir: str | Path, image_id: str) -> dict[str, Path]:
    out_dir = Path(out_dir)
    return {k: out_dir / f"{image_id}{suf}" for k, suf in S.SUFFIXES.items()}


# ------------------------------------------------------------------------------------ overview
CLASS_COLOURS = {"complete": "#2ca02c", "partial": "#1f77b4", "merged_suspect": "#d62728", "fragment": "#aaaaaa",
                 "verified": "#9467bd"}


def overview_png(path: str | Path, ce: np.ndarray, strands: pd.DataFrame, samples: pd.DataFrame,
                 nuclei: pd.DataFrame, nuclei_run: pd.DataFrame, spacing, image_id: str, stride: int = 2) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mip = np.asarray(ce).max(axis=0)[::stride, ::stride].astype(np.float32)
    lo, hi = np.percentile(mip, [1, 99.7])
    dy, dx = float(spacing[1]) * stride, float(spacing[2]) * stride
    h, w = mip.shape
    fig = plt.figure(figsize=(10, 10 * h / max(w, 1)), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(mip, cmap="gray", vmin=lo, vmax=hi, extent=[0, w * dx, h * dy, 0], interpolation="nearest")
    if len(samples):
        cls = strands.set_index("strand_id")["strand_class"].to_dict() if len(strands) else {}
        smp = samples[samples["in_strand"].astype(bool)]
        for sid, g in smp.groupby("strand_id"):
            ax.plot(g["x_um"], g["y_um"], "-", lw=0.6, color=CLASS_COLOURS.get(cls.get(sid, ""), "#ff7f0e"))
    if len(nuclei):
        bad = nuclei[nuclei["status"] == "failed"]["nucleus_id"]
        cen = nuclei_run.set_index("nucleus_id").reindex(bad)
        if len(cen):
            ax.plot(cen["centroid_x_um"], cen["centroid_y_um"], "o", mfc="none", mec="red", ms=10)
    counts = strands["strand_class"].value_counts().to_dict() if len(strands) else {}
    st = nuclei["status"].value_counts().to_dict() if len(nuclei) else {}
    txt = (f"{image_id}\nnuclei: " + ", ".join(f"{k} {v}" for k, v in sorted(st.items())) +
           "\nstrands: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    ax.text(0.01, 0.01, txt, transform=ax.transAxes, fontsize=8, color="white", va="bottom",
            bbox={"facecolor": "black", "alpha": 0.6, "pad": 3})
    for k, c in CLASS_COLOURS.items():
        ax.plot([], [], "-", color=c, label=k)
    ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    ax.set_xlim(0, w * dx)
    ax.set_ylim(h * dy, 0)
    ax.set_axis_off()
    path = Path(path)
    tmp = _tmp(path)
    fig.savefig(long_path(tmp), dpi=100, format="png")
    plt.close(fig)
    _replace(tmp, path)


# ------------------------------------------------------------------------------------ collect
def _records(out_root: Path) -> list[tuple[Path, dict]]:
    root_lp = long_path(out_root)
    hits = []
    for dirpath, _dirs, names in os.walk(root_lp):
        rel = os.path.relpath(dirpath, root_lp)
        base = out_root if rel == "." else out_root / rel
        for n in names:
            if n.endswith(S.SUFFIXES["record"]) and not n.startswith("batch_"):
                rec = read_record(base / n)
                if rec and rec.get("format") == "germquant.sc_straighten":
                    hits.append((base, rec))
    return sorted(hits, key=lambda t: (str(t[0]), str(t[1].get("image_id"))))


def collect_tables(out_root: str | Path) -> dict | None:
    """Stack the sc_straighten tables of every finished record under ``out_root`` into ``batch_sc_*.csv`` and
    write ``batch_sc_straighten_summary.csv``. Returns {name: rows}, or None (and writes nothing) when the tree
    holds no sc_straighten record."""
    out_root = Path(out_root)
    recs = _records(out_root)
    if not recs:
        return None
    counts = {}
    done = [(d, r) for d, r in recs if r.get("done")]
    for t in S.COLLECTED:
        frames = []
        for d, r in done:
            f = d / f"{r['image_id']}{S.SUFFIXES[t]}"
            if f.is_file():
                try:
                    frames.append(pd.read_csv(long_path(f), low_memory=False))
                except Exception as e:  # noqa: BLE001
                    log.warning("collect: skipping %s (%s)", f, e)
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=S.TABLES[t])
        df.to_csv(long_path(out_root / S.BATCH_NAMES[t]), index=False, float_format="%.6g", lineterminator="\n")
        counts[f"sc_{t}"] = len(df)
    rows = []
    for d, r in recs:
        c = r.get("counts", {})
        rows.append({"image_id": r.get("image_id"), "status": r.get("status"),
                     "sc_straighten_version": r.get("version"), "sc_straighten_hash": r.get("sub_hash"),
                     "params_hash": r.get("params_hash"), "n_nuclei": c.get("nuclei"),
                     "n_traced": c.get("status_traced", 0), "n_no_sc": c.get("status_no_sc", 0),
                     "n_failed": c.get("status_failed", 0), "n_skipped": c.get("status_skipped", 0),
                     "n_strands": c.get("strands"), "n_complete": c.get("class_complete", 0),
                     "n_partial": c.get("class_partial", 0), "n_merged_suspect": c.get("class_merged_suspect", 0),
                     "n_fragment": c.get("class_fragment", 0), "n_verified": c.get("class_verified", 0),
                     "flags": ";".join(r.get("flags", [])), "elapsed_s": r.get("elapsed_s"), "results_dir": str(d)})
    pd.DataFrame(rows, columns=S.SUMMARY).to_csv(long_path(out_root / S.BATCH_NAMES["summary"]), index=False,
                                                 lineterminator="\n")
    counts["sc_straighten_summary"] = len(rows)
    return counts
