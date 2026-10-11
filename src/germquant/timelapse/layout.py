"""What a finished series folder looks like from the outside, read-only (docs/TIMELAPSE.md, "Tracking").

Every tracking command (``germquant track-export``, ``track-link``, ``track-import``, ``track-strips``,
``track-validate`` and scripts/pick_validation_nuclei.py) starts from the folder ``process_series`` wrote:
``<stem>__series.json`` next to one ``t<NNN>/`` folder per timepoint, each holding the fixed-image file set
under the id ``<stem>_t<NNN>``. ``SeriesLayout`` finds the stem, the timepoints that are actually on disk
(a nuclei table and a label image), their recorded acquisition times, the voxel size and, when the
registration package wrote it, ``<stem>__register_shifts.csv``. Nothing here writes to the series folder.

Registered coordinates (the shared contract of the time-lapse packages): the raw centroid of a nucleus at
timepoint k minus the cumulative shift of that timepoint, ``reg = raw - shift_k``, where a positive shift
means the sample moved by that much between t0 and k. Nothing is ever resampled: the per-timepoint tables
and label images stay in their own raw frame, and every tracking step works on coordinates.
"""
from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .. import schema
from ..fsutil import long_path
from ..provenance import read_done_marker
from .series import SERIES_RECORD_SUFFIX, read_series_record

log = logging.getLogger(__name__)

SHIFTS_TABLE = "register_shifts"
# the cumulative-shift columns of <stem>__register_shifts.csv: the shared-contract names first, then the
# names schema.REGISTER_SHIFTS declared before the contract was fixed (a table written with either set is read)
_SHIFT_COLUMNS = (("dz_um", "dy_um", "dx_um"), ("cum_shift_z_um", "cum_shift_y_um", "cum_shift_x_um"))
_REF_RE = re.compile(r"__ref_stride(\d+)\.npy$")
_T_DIR_RE = re.compile(r"^t(\d{3,})$")
_ID_RE = re.compile(r"^(?P<stem>.+)_t(?P<t>\d{3,})$")

# nucleus measurements joined onto the track points when the nuclei table carries them (schema.NUCLEI
# names plus the per-role intensity extras the measure stage writes)
JOIN_COLUMNS = (
    "volume_um3", "n_voxels", "in_germline", "axis_position_norm", "axis_position_um", "zone", "n_spots",
    "sc_total_length_um", "sc_n_fragments_lb", "sc_fragmentation_index", "sc_expected_n_tracks",
)
JOIN_INTENSITY_ROLES = ("dna", "central_element")


class SeriesLayoutError(RuntimeError):
    """The folder is not a finished series folder this package can read."""


@dataclass
class SeriesLayout:
    out_dir: Path
    stem: str
    record: dict | None
    timepoints: list[int]                      # the timepoints with a nuclei table and a label image on disk
    file_path: str | None = None               # the series file (from the record), for `--image file`
    period_s: float | None = None
    n_timepoints: int | None = None
    t_seconds: dict[int, float | None] = field(default_factory=dict)
    missing: list[int] = field(default_factory=list)   # timepoints the record lists as done that are not on disk

    # ---- paths -------------------------------------------------------------------------------------

    def sub(self, t: int) -> Path:
        return self.out_dir / f"t{int(t):03d}"

    def image_id(self, t: int) -> str:
        return f"{self.stem}_t{int(t):03d}"

    def nuclei_path(self, t: int) -> Path:
        return self.sub(t) / f"{self.image_id(t)}__nuclei.csv"

    def labels_path(self, t: int) -> Path:
        return self.sub(t) / f"{self.image_id(t)}__nuclei_labels.tif"

    def summary_path(self, t: int) -> Path:
        return self.sub(t) / f"{self.image_id(t)}__image_summary.csv"

    def ref_path(self, t: int) -> tuple[Path, int] | None:
        """The downsampled copy of the segmentation channel the loop keeps for the registration and tracking
        packages, ``t<NNN>/<id>__ref_stride<k>.npy`` (uint16, Z x Y/k x X/k), and its stride k; None when the
        series was run without it."""
        sub = self.sub(t)
        try:
            names = sorted(p.name for p in sub.iterdir())
        except OSError:
            return None
        prefix = f"{self.image_id(t)}__ref_stride"
        for name in names:
            if name.startswith(prefix):
                m = _REF_RE.search(name)
                if m:
                    return sub / name, int(m.group(1))
        return None

    def series_table_path(self, name: str) -> Path:
        return self.out_dir / f"{self.stem}__{name}.csv"

    # ---- readers -----------------------------------------------------------------------------------

    def read_nuclei(self, t: int) -> pd.DataFrame:
        df = pd.read_csv(long_path(self.nuclei_path(t)), low_memory=False)
        if "nucleus_id" in df.columns:
            df["nucleus_id"] = pd.to_numeric(df["nucleus_id"], errors="coerce")
            df = df[df["nucleus_id"].notna()].copy()
            df["nucleus_id"] = df["nucleus_id"].astype(int)
        return df

    def read_labels(self, t: int) -> np.ndarray:
        import tifffile

        lab = np.asarray(tifffile.imread(long_path(self.labels_path(t))))
        if lab.ndim == 2:
            lab = lab[None]
        if lab.ndim != 3:
            raise SeriesLayoutError(f"{self.labels_path(t).name}: expected a (Z, Y, X) label image, got shape {lab.shape}")
        return lab

    def labels_shape(self, t: int) -> tuple[int, int, int]:
        """The (Z, Y, X) shape of a timepoint's label image without reading its pixels."""
        import tifffile

        with tifffile.TiffFile(long_path(self.labels_path(t))) as tf:
            shape = tuple(int(s) for s in tf.series[0].shape)
        if len(shape) == 2:
            shape = (1, *shape)
        return shape  # type: ignore[return-value]

    def read_summary(self, t: int) -> pd.Series | None:
        p = self.summary_path(t)
        try:
            s = pd.read_csv(long_path(p), low_memory=False)
        except (OSError, ValueError):
            return None
        return s.iloc[0] if len(s) else None

    def spacing(self) -> tuple[float, float, float]:
        """(dz, dy, dx) in microns from the voxel columns of the first nuclei table that has a row."""
        for t in self.timepoints:
            df = self.read_nuclei(t)
            if len(df) and {"voxel_dz_um", "voxel_dy_um", "voxel_dx_um"} <= set(df.columns):
                r = df.iloc[0]
                sp = (float(r["voxel_dz_um"]), float(r["voxel_dy_um"]), float(r["voxel_dx_um"]))
                if all(np.isfinite(sp)) and all(v > 0 for v in sp):
                    return sp
        raise SeriesLayoutError(f"{self.stem}: no nuclei table carries a voxel size (voxel_dz_um, voxel_dy_um, "
                                f"voxel_dx_um); every timepoint is empty?")

    def run_xy_stride(self) -> int:
        """The xy stride the series was run at (from the series record; 1 when the record has none). The label
        images and the nuclei tables are on that grid, the series file is not."""
        try:
            return max(1, int((self.record or {}).get("xy_stride") or 1))
        except (TypeError, ValueError):
            return 1

    def run_z_range(self) -> tuple[int, int] | None:
        """The z range the series was run on (from the series record), or None for every plane."""
        zr = (self.record or {}).get("z_range")
        return (int(zr[0]), int(zr[1])) if zr else None

    def channel_index(self, role: str) -> int | None:
        """The channel index a role resolved to, from the channel_map_json column of the first nuclei table."""
        for t in self.timepoints:
            df = self.read_nuclei(t)
            if len(df) and "channel_map_json" in df.columns:
                try:
                    m = json.loads(str(df["channel_map_json"].iloc[0]))
                except ValueError:
                    continue
                idx = m.get(role)
                return int(idx) if idx is not None else None
        return None

    def t_seconds_of(self, t: int) -> float | None:
        v = self.t_seconds.get(int(t))
        if v is None and self.period_s:
            return float(t) * float(self.period_s)
        return None if v is None else float(v)

    def minutes_of(self, t: int) -> float:
        """Minutes since the first timepoint on disk (the recorded stamps when the file has them, else
        t * period_s, else the index itself: one unit per timepoint)."""
        t0 = self.timepoints[0] if self.timepoints else 0
        a, b = self.t_seconds_of(t), self.t_seconds_of(t0)
        if a is not None and b is not None:
            return (a - b) / 60.0
        return float(int(t) - int(t0))

    def series_meta(self, run_timestamp: str | None = None) -> dict:
        """The SERIES_META block for write_series_table: the design columns from the first nuclei table,
        the counts from the record, the provenance of this process."""
        from .. import provenance

        meta = {k: None for k in schema.SERIES_META}
        meta.update(series_id=self.stem, file_path=self.file_path,
                    n_timepoints=self.n_timepoints, period_s=self.period_s)
        for t in self.timepoints:
            df = self.read_nuclei(t)
            if len(df):
                r = df.iloc[0]
                for k in ("genotype", "sex", "germ_cell", "treatment", "replicate", "acquisition_date"):
                    if k in df.columns:
                        meta[k] = r[k]
                if meta["file_path"] is None and "file_path" in df.columns:
                    meta["file_path"] = r["file_path"]
                break
        rec = self.record or {}
        meta["pipeline_version"] = rec.get("pipeline_version") or provenance._safe_version("germquant")
        meta["git_sha"] = provenance.git_sha()
        meta["config_hash"] = rec.get("config_hash")
        meta["run_timestamp"] = run_timestamp or provenance.run_timestamp()
        return meta


def find_series(out_dir: str | Path, stem: str | None = None) -> SeriesLayout:
    """Describe the series folder ``out_dir``. The stem is the given one, else the single ``<stem>__series.json``
    in the folder, else the one stem the ``t<NNN>/<stem>_t<NNN>__nuclei.csv`` files agree on. The timepoints
    are those with both a nuclei table and a label image; when a record exists, only the timepoints it lists
    as done are used and the others are reported in ``missing`` (a failed or interrupted timepoint is a gap in
    the series, not a frame of zeros)."""
    out_dir = Path(out_dir)
    if not out_dir.is_dir():
        raise SeriesLayoutError(f"{out_dir}: not a folder")
    names = sorted(p.name for p in out_dir.iterdir())
    if stem is None:
        records = [n[: -len(SERIES_RECORD_SUFFIX)] for n in names if n.endswith(SERIES_RECORD_SUFFIX)]
        if len(records) == 1:
            stem = records[0]
        elif len(records) > 1:
            raise SeriesLayoutError(f"{out_dir}: several series records ({', '.join(records)}); pass the stem")
        else:
            stems = set()
            for n in names:
                if _T_DIR_RE.match(n) and (out_dir / n).is_dir():
                    for f in (out_dir / n).glob("*__nuclei.csv"):
                        m = _ID_RE.match(f.name.split("__")[0])
                        if m:
                            stems.add(m.group("stem"))
            if len(stems) != 1:
                raise SeriesLayoutError(f"{out_dir}: no <stem>__series.json and {len(stems)} candidate stems in the "
                                        f"t<NNN> folders ({', '.join(sorted(stems)) or 'none'}); pass the stem")
            stem = stems.pop()
    record = read_series_record(out_dir, stem)
    layout = SeriesLayout(out_dir=out_dir, stem=stem, record=record, timepoints=[])
    on_disk = []
    for n in names:
        m = _T_DIR_RE.match(n)
        if not m or not (out_dir / n).is_dir():
            continue
        k = int(m.group(1))
        if layout.nuclei_path(k).is_file() and layout.labels_path(k).is_file():
            on_disk.append(k)
    on_disk.sort()
    if record:
        done = sorted(int(k) for k in record.get("done", []))
        layout.timepoints = [k for k in done if k in set(on_disk)] if done else on_disk
        layout.missing = [k for k in done if k not in set(on_disk)]
        extra = [k for k in on_disk if k not in set(done)]
        if extra:
            log.warning("%s: timepoint folder(s) %s are on disk but the series record does not list them as done; "
                        "they are left out", stem, extra)
        layout.file_path = record.get("file_path")
        layout.period_s = record.get("period_s")
        layout.n_timepoints = record.get("n_timepoints")
        stamps = record.get("t_seconds")
        for e in record.get("timepoints", []):
            if "t_index" in e:
                layout.t_seconds[int(e["t_index"])] = e.get("t_seconds")
        if isinstance(stamps, list):
            for k, v in enumerate(stamps):
                layout.t_seconds.setdefault(k, v)
    else:
        layout.timepoints = on_disk
        layout.n_timepoints = (max(on_disk) + 1) if on_disk else None
        for k in on_disk:                     # the completion markers carry the stamp when the file has one
            marker = read_done_marker(layout.sub(k), layout.image_id(k)) or {}
            layout.t_seconds[k] = marker.get("t_seconds")
    if layout.missing:
        log.warning("%s: timepoint(s) %s are in the series record but their tables are not on disk", stem, layout.missing)
    if not layout.timepoints:
        raise SeriesLayoutError(f"{out_dir}: no timepoint of {stem} has both a nuclei table and a label image")
    return layout


def load_shifts(layout: SeriesLayout) -> tuple[pd.DataFrame, str]:
    """The cumulative shift to t0 per timepoint (``t_index``, ``dz_um``, ``dy_um``, ``dx_um``; microns, raw
    frame) from ``<stem>__register_shifts.csv``, and where it came from: ``"register_shifts"`` or ``"none"``
    (no table: zero shifts, the series is treated as unregistered). A timepoint the table lacks, or a NaN
    shift, is a zero shift with a warning; a positive shift means the sample moved by that much."""
    zero = pd.DataFrame({"t_index": layout.timepoints, "dz_um": 0.0, "dy_um": 0.0, "dx_um": 0.0})
    p = layout.series_table_path(SHIFTS_TABLE)
    if not p.is_file():
        log.info("%s: no %s; coordinates are taken as already aligned (zero shifts)", layout.stem, p.name)
        return zero, "none"
    df = pd.read_csv(long_path(p), low_memory=False)
    if "t_index" not in df.columns:
        raise SeriesLayoutError(f"{p.name}: no t_index column")
    cols = next((c for c in _SHIFT_COLUMNS if set(c) <= set(df.columns)), None)
    if cols is None:
        raise SeriesLayoutError(f"{p.name}: no cumulative shift columns (expected {', '.join(_SHIFT_COLUMNS[0])})")
    df = df[["t_index", *cols]].copy()
    df.columns = ["t_index", "dz_um", "dy_um", "dx_um"]
    df["t_index"] = pd.to_numeric(df["t_index"], errors="coerce")
    df = df[df["t_index"].notna()].drop_duplicates("t_index", keep="last")
    df["t_index"] = df["t_index"].astype(int)
    out = zero.merge(df, on="t_index", how="left", suffixes=("_zero", ""))
    for c in ("dz_um", "dy_um", "dx_um"):
        v = pd.to_numeric(out[c], errors="coerce")
        bad = v.isna()
        if bad.any():
            log.warning("%s: %d timepoint(s) without a %s value; zero shift used for %s", p.name, int(bad.sum()), c,
                        out.loc[bad, "t_index"].tolist())
        out[c] = v.fillna(0.0).astype(float)
    return out[["t_index", "dz_um", "dy_um", "dx_um"]].reset_index(drop=True), "register_shifts"


def shift_lookup(shifts: pd.DataFrame) -> dict[int, np.ndarray]:
    """``{t_index: array([dz, dy, dx]) um}`` from load_shifts' table."""
    return {int(r.t_index): np.array([r.dz_um, r.dy_um, r.dx_um], dtype=float) for r in shifts.itertuples(index=False)}


def voxel_shifts(shifts: pd.DataFrame, spacing) -> dict[int, np.ndarray]:
    """The cumulative shift per timepoint rounded to whole voxels (``round(shift / spacing)``), the amount a
    label image is displaced by when it is placed on the registered grid."""
    sp = np.asarray(spacing, dtype=float)
    return {t: np.rint(v / sp).astype(int) for t, v in shift_lookup(shifts).items()}


def object_table(layout: SeriesLayout, shifts: pd.DataFrame, *, germline_only: bool = True) -> pd.DataFrame:
    """One row per nucleus per timepoint from the per-timepoint nuclei tables: ``t_index``, ``t_seconds``,
    ``image_id``, ``nucleus_id``, the raw centroid (``z_um``, ``y_um``, ``x_um``), the registered centroid
    (``zr_um``, ``yr_um``, ``xr_um`` = raw minus the timepoint's cumulative shift) and the measurements the
    track points carry (JOIN_COLUMNS and the dna / central_element mean intensities when present). With
    ``germline_only`` the rows with in_germline False are dropped (a table without the column keeps every
    nucleus, as the germline stage being off means every nucleus is the analysis set)."""
    sh = shift_lookup(shifts)
    frames = []
    for t in layout.timepoints:
        df = layout.read_nuclei(t)
        if not len(df):
            continue
        if germline_only and "in_germline" in df.columns:
            keep = df["in_germline"].map(_truthy)
            df = df[keep].copy()
        out = pd.DataFrame({
            "t_index": int(t), "t_seconds": layout.t_seconds_of(t), "image_id": layout.image_id(t),
            "nucleus_id": df["nucleus_id"].astype(int).to_numpy(),
            "z_um": pd.to_numeric(df["centroid_z_um"], errors="coerce").to_numpy(float),
            "y_um": pd.to_numeric(df["centroid_y_um"], errors="coerce").to_numpy(float),
            "x_um": pd.to_numeric(df["centroid_x_um"], errors="coerce").to_numpy(float),
        })
        s = sh.get(int(t), np.zeros(3))
        out["zr_um"] = out["z_um"] - s[0]
        out["yr_um"] = out["y_um"] - s[1]
        out["xr_um"] = out["x_um"] - s[2]
        for c in JOIN_COLUMNS:
            if c in df.columns:
                out[c] = df[c].to_numpy()
        for role in JOIN_INTENSITY_ROLES:
            c = f"{role}_mean_intensity"
            if c in df.columns:
                out[c] = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
        frames.append(out)
    if not frames:
        cols = ["t_index", "t_seconds", "image_id", "nucleus_id", "z_um", "y_um", "x_um", "zr_um", "yr_um", "xr_um"]
        return pd.DataFrame(columns=cols)
    return pd.concat(frames, ignore_index=True)


def _truthy(v) -> bool:
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes")
    try:
        return bool(v) and not pd.isna(v)
    except (TypeError, ValueError):
        return False
