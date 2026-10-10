"""The per-timepoint loop around process_image (docs/DESIGN_timelapse.md section 2, docs/TIMELAPSE.md).

``process_series`` runs every enabled stage on each requested timepoint of one time-lapse file, exactly as a
fixed image is run, into ``out_dir/t<NNN>/``: each of those folders holds the fixed-image file set under
the image id ``<stem>_t<NNN>`` (process_image with ``t``). Next to them it keeps two series-level records:

* ``<stem>__series.json``: the timepoints requested, done and failed (with their flags), period_s and
  t_seconds from the header, the elapsed time per timepoint, config_hash and pipeline version. It is
  rewritten after every timepoint, so an interrupted series leaves an honest record, and ``resume=True``
  skips the timepoints whose completion marker matches (the rule ``germquant batch --resume`` uses) and
  reprocesses the rest.
* ``<stem>__series_summary.csv`` (and .parquet): one row per requested timepoint from that timepoint's
  image_summary (schema.SERIES_SUMMARY plus the signal_qc background and SNR per channel role).

A timepoint whose run raises is recorded (an ``EXCEPTION:<type>: <text>`` flag, no completion marker) and
the loop goes on, as a batch does. Nothing here changes a per-image table or the fixed-image path.
"""
from __future__ import annotations

import json
import logging
import os
import time
from pathlib import Path

import pandas as pd

from .. import batch as B
from .. import pipeline, provenance, schema
from ..fsutil import long_path
from ..io import parse_sample
from ..stages import stage_hashes as stage_hashes_fn

log = logging.getLogger(__name__)

SERIES_RECORD_SUFFIX = "__series.json"
SUMMARY_TABLE = "series_summary"


def parse_timepoints(spec, n_timepoints: int) -> list[int]:
    """The timepoint indices a ``--t`` value names, in increasing order: ``"all"`` (every timepoint of the
    file), one index (``7`` or ``"7"``), an inclusive range ``"a-b"`` (``"0-24"``), or an iterable of
    indices. ValueError when an index is outside ``0 .. n_timepoints - 1``, when a range runs backwards, or
    when the text is none of these."""
    n = int(n_timepoints)
    if n < 1:
        raise ValueError(f"the file has no timepoints (n_timepoints={n})")

    def _check(k) -> int:
        try:
            k = int(k)
        except (TypeError, ValueError) as e:
            raise ValueError(f"timepoint {k!r} is not an integer") from e
        if not 0 <= k < n:
            raise ValueError(f"timepoint {k} is out of range: the file has {n} timepoint(s), 0..{n - 1}")
        return k

    if isinstance(spec, bool):
        raise ValueError(f"timepoints must be 'all', an index, 'a-b' or a list of indices, got {spec!r}")
    if isinstance(spec, int):
        return [_check(spec)]
    if isinstance(spec, str):
        text = spec.strip().lower()
        if text == "all":
            return list(range(n))
        if "-" in text[1:]:                      # "a-b" (a leading minus is not a range)
            a, b = text.split("-", 1)
            lo, hi = _check(a), _check(b)
            if hi < lo:
                raise ValueError(f"timepoint range {spec!r} runs backwards")
            return list(range(lo, hi + 1))
        return [_check(text)]
    try:
        items = list(spec)
    except TypeError as e:
        raise ValueError(f"timepoints must be 'all', an index, 'a-b' or a list of indices, got {spec!r}") from e
    out = sorted({_check(k) for k in items})
    if not out:
        raise ValueError("no timepoints requested")
    return out


def read_series_record(out_dir: str | Path, stem: str) -> dict | None:
    """The ``<stem>__series.json`` record under out_dir, or None when there is none (or it is unreadable)."""
    p = Path(out_dir) / f"{stem}{SERIES_RECORD_SUFFIX}"
    try:
        with open(long_path(p), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def process_series(
    path: str | Path,
    cfg,
    out_dir: str | Path,
    *,
    timepoints="all",
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    resume: bool = False,
) -> dict:
    """Run process_image on each requested timepoint of ``path`` (``timepoints``: "all", an index, "a-b" or a
    list; see parse_timepoints) into ``out_dir/t<NNN>/`` and write the series record and the series_summary
    table next to those folders (module docstring). ``xy_stride`` and ``z_range`` are handed to every
    timepoint. With ``resume`` the timepoints whose completion marker matches this config, geometry and
    segment stage are skipped and their summary rows are read back from disk.

    Returns ``series_id``, ``out_dir``, ``n_timepoints`` (of the file), ``requested``, ``done`` and ``failed``
    (t_index lists; ``failed`` as ``{"t_index", "flags"}`` dicts), ``record`` (the series.json content),
    ``record_path``, ``summary`` (the series_summary DataFrame) and ``summary_path``."""
    path = Path(path)
    out_dir = Path(out_dir)
    xy_stride = int(xy_stride) if xy_stride else 1
    zr = (int(z_range[0]), int(z_range[1])) if z_range else None
    os.makedirs(long_path(out_dir), exist_ok=True)
    started = provenance.run_timestamp()

    # the header through the pipeline's reader names (the seam the synthetic test fixtures patch)
    meta = pipeline.read_input_metadata(path)
    n_t = int(meta.get("n_timepoints") or (meta.get("sizes") or {}).get("T", 1) or 1)
    period_s = meta.get("period_s")
    t_seconds = meta.get("t_seconds")
    requested = parse_timepoints(timepoints, n_t)

    sample = parse_sample(path, cfg.get("metadata.filename_regex"), cfg.get("metadata.defaults"))
    stem = sample["image_id"]                     # the raw stem: process_image suffixes it per timepoint

    # series-level provenance, as `germquant batch` keeps it: one manifest at the root and the segment
    # sub-hash that lets --resume notice a retrained model or a Cellpose upgrade
    model_sha = (provenance.file_sha256(pipeline.resolve_model_path(cfg))
                 if str(cfg.get("segmentation.nuclei.method", "auto")) != "classical" else "unused")
    segment_hash = stage_hashes_fn(cfg, {}, provenance.tool_versions(), model_sha)["segment"]
    enabled = B.enabled_stage_names(cfg)
    prov = provenance.write_manifest(
        out_dir, config_hash=cfg.hash, config=cfg.as_dict(),
        extra={"series_id": stem, "input": str(path), "n_timepoints": n_t, "timepoints_requested": requested,
               "period_s": period_s, "xy_stride": xy_stride, "z_range": list(zr) if zr else None,
               "model_sha256": model_sha, "enabled_stages": enabled, "resume": bool(resume)})

    prev = read_series_record(out_dir, stem) if resume else None
    prev_entries = {int(e["t_index"]): e for e in (prev or {}).get("timepoints", []) if "t_index" in e}

    record = {
        "series_id": stem, "file_path": str(path),
        "config_hash": cfg.hash, "pipeline_version": prov["pipeline_version"], "git_sha": prov["git_sha"],
        "n_timepoints": n_t, "period_s": period_s, "t_seconds": t_seconds,
        "xy_stride": xy_stride, "z_range": list(zr) if zr else None, "enabled_stages": enabled,
        "resume": bool(resume), "requested": requested, "done": [], "failed": [], "elapsed_s": {},
        "timepoints": [], "started": started, "finished": None, "complete": False,
    }
    record_path = out_dir / f"{stem}{SERIES_RECORD_SUFFIX}"
    rows: list[dict] = []

    for i, k in enumerate(requested, 1):
        sub = out_dir / f"t{k:03d}"
        iid = f"{stem}_t{k:03d}"
        ts = pipeline.t_seconds_of(meta, k)
        entry = {"t_index": k, "t_seconds": ts, "image_id": iid, "out_dir": sub.name, "status": "",
                 "resumed": False, "qc_pass": None, "qc_flags": [], "failed_stages": [], "elapsed_s": None}
        if resume:
            done, why = B.is_done(sub, iid, cfg, xy_stride=xy_stride, z_range=zr, segment_hash=segment_hash, t=k)
            if done:
                marker = provenance.read_done_marker(sub, iid) or {}
                row = _summary_row(sub, iid, k, ts)
                entry.update(status="done", resumed=True, qc_pass=bool(row["qc_pass"]),
                             qc_flags=[f for f in str(row["qc_flags"]).split(";") if f],
                             failed_stages=list(marker.get("failed_stages", [])),
                             elapsed_s=prev_entries.get(k, {}).get("elapsed_s"))
                row["elapsed_s"] = entry["elapsed_s"]
                rows.append(row)
                log.info("[%d/%d] %s already done (completion marker matches); skipping", i, len(requested), iid)
                _append(record, entry, record_path)
                continue
            if why != "no completion marker":
                log.info("[%d/%d] %s reprocessing: %s", i, len(requested), iid, why)
        log.info("[%d/%d] %s", i, len(requested), iid)
        t0 = time.perf_counter()
        try:
            res = pipeline.process_image(path, cfg, sub, xy_stride=xy_stride, z_range=zr, prov=None, t=k)
            elapsed = round(time.perf_counter() - t0, 3)
            row = _summary_row(sub, res["image_id"], k, ts)
            row["elapsed_s"] = elapsed
            entry.update(status="done", qc_pass=bool(res["qc_pass"]), qc_flags=list(res["qc_flags"]),
                         failed_stages=[n for n, oc in res["stages"].items() if oc.get("status") == "failed"],
                         elapsed_s=elapsed)
        except Exception as e:  # noqa: BLE001 - one timepoint failing must not stop the series
            elapsed = round(time.perf_counter() - t0, 3)
            log.exception("%s FAILED", iid)
            flag = f"EXCEPTION:{type(e).__name__}: {e}"
            entry.update(status="failed", qc_pass=False, qc_flags=[flag], elapsed_s=elapsed)
            row = _failed_row(iid, k, ts, flag, elapsed)
        rows.append(row)
        _append(record, entry, record_path)

    record["finished"] = provenance.run_timestamp()
    record["complete"] = True
    _write_json(record_path, record)

    series_meta = {
        "series_id": stem, "file_path": str(path),
        "genotype": sample["genotype"], "sex": sample["sex"], "germ_cell": sample["germ_cell"],
        "treatment": sample["treatment"], "replicate": sample["replicate"],
        "acquisition_date": sample["acquisition_date"],
        "n_timepoints": n_t, "period_s": period_s,
        "pipeline_version": prov["pipeline_version"], "git_sha": prov["git_sha"],
        "config_hash": cfg.hash, "run_timestamp": prov["run_timestamp"],
    }
    summary = _summary_table(rows)
    summary_path = write_series_table(summary, SUMMARY_TABLE, out_dir, stem, series_meta)

    n_done, n_failed = len(record["done"]), len(record["failed"])
    log.info("%s: %d timepoint(s) requested, %d done, %d failed -> %s", stem, len(requested), n_done, n_failed, out_dir)
    return {"series_id": stem, "out_dir": str(out_dir), "n_timepoints": n_t, "requested": requested,
            "done": list(record["done"]), "failed": [dict(f) for f in record["failed"]],
            "record": record, "record_path": str(record_path),
            "summary": summary, "summary_path": str(summary_path)}


def write_series_table(df: pd.DataFrame, name: str, out_dir: str | Path, stem: str, series_meta: dict) -> Path:
    """Write one series table (schema.SERIES_TABLES) as ``<stem>__<name>.csv`` and ``.parquet`` under
    out_dir: the declared columns first (NA where absent), extras after them, then the SERIES_META block on
    every row (pipeline._write_tables' rule for the per-image tables). The parquet copy is best effort (a
    warning when pyarrow is missing); the CSV path is returned."""
    out_dir = Path(out_dir)
    declared = list(schema.SERIES_TABLES[name])
    df = df.copy()
    for col in declared:
        if col not in df.columns:
            df[col] = pd.NA
    extras = [c for c in df.columns if c not in declared and c not in schema.SERIES_META]
    df = df[declared + extras]
    for k in schema.SERIES_META:
        df[k] = series_meta.get(k)
    base = out_dir / f"{stem}__{name}"
    csv = base.with_suffix(".csv")
    df.to_csv(long_path(csv), index=False)
    try:
        df.to_parquet(long_path(base.with_suffix(".parquet")), index=False)
    except Exception as e:  # noqa: BLE001 - pyarrow missing or a dtype it refuses: the CSV is the record
        log.warning("parquet write failed for %s (%s); CSV written.", base.name, e)
    return csv


# ---- helpers ------------------------------------------------------------------------------------------

def _append(record: dict, entry: dict, record_path: Path) -> None:
    """Add a timepoint's entry to the record and rewrite the series.json (after every timepoint, so an
    interrupted series leaves an honest record of what is on disk)."""
    record["timepoints"].append(entry)
    k = int(entry["t_index"])
    if entry["status"] == "done":
        record["done"] = sorted(set(record["done"]) | {k})
    elif entry["status"] == "failed":
        record["failed"].append({"t_index": k, "flags": list(entry["qc_flags"])})
    record["elapsed_s"][str(k)] = entry["elapsed_s"]
    _write_json(record_path, record)


def _write_json(path: Path, obj: dict) -> None:
    with open(long_path(path), "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, default=str)


def _roles_present(row: pd.Series) -> list[str]:
    """The channel roles the timepoint resolved, in channel-map order, from the image_summary row's
    channel_map_json; falls back to the roles whose sig_<role>_bg is a number."""
    cm = row.get("channel_map_json")
    if isinstance(cm, str) and cm:
        try:
            m = json.loads(cm)
            return [r for r, idx in m.items() if idx is not None]
        except ValueError:
            pass
    roles = list(schema.SIGNAL_QC_ROLES)
    roles += sorted({c[len("sig_"):-len("_bg")] for c in row.index
                     if c.startswith("sig_") and c.endswith("_bg") and c[len("sig_"):-len("_bg")] not in roles})
    return [r for r in roles if pd.notna(row.get(f"sig_{r}_bg"))]


def _summary_row(sub: Path, image_id: str, t_index: int, t_seconds) -> dict:
    """The series_summary row of a finished timepoint, read back from its ``<image_id>__image_summary.csv``
    (the same source for a timepoint processed now and for one a resume skipped)."""
    f = sub / f"{image_id}__image_summary.csv"
    s = pd.read_csv(long_path(f), low_memory=False, float_precision="round_trip")
    r = s.iloc[0] if not s.empty else pd.Series(dtype=object)
    flags = r.get("qc_flags", "")
    row = {
        "t_index": int(t_index), "t_seconds": t_seconds, "image_id": image_id,
        "n_nuclei": B._int(r.get("n_nuclei"), 0), "n_germline_nuclei": B._int(r.get("n_germline_nuclei"), 0),
        "mean_spots": float(r["mean_spots"]) if pd.notna(r.get("mean_spots")) else float("nan"),
        "qc_pass": B._bool(r.get("qc_pass", False)),
        "qc_flags": str(flags) if pd.notna(flags) else "",
        "elapsed_s": None,
    }
    for role in _roles_present(r):
        for metric in ("bg", "snr"):
            col = f"sig_{role}_{metric}"
            row[col] = float(r[col]) if pd.notna(r.get(col)) else float("nan")
    return row


def _failed_row(image_id: str, t_index: int, t_seconds, flag: str, elapsed: float) -> dict:
    return {"t_index": int(t_index), "t_seconds": t_seconds, "image_id": image_id,
            "n_nuclei": pd.NA, "n_germline_nuclei": pd.NA, "mean_spots": float("nan"),
            "qc_pass": False, "qc_flags": flag, "elapsed_s": elapsed}


def _summary_table(rows: list[dict]) -> pd.DataFrame:
    """One row per requested timepoint, in t_index order: the declared columns, then the sig_<role>_bg /
    sig_<role>_snr columns in the order they were first seen (the channel map's role order)."""
    cols = list(schema.SERIES_SUMMARY)
    for r in rows:
        cols += [c for c in r if c not in cols]
    df = pd.DataFrame(rows, columns=cols)
    if len(df):
        df = df.sort_values("t_index", kind="stable").reset_index(drop=True)
    df["t_seconds"] = pd.to_numeric(df["t_seconds"], errors="coerce")      # NaN when the header has no times
    for c in ("n_nuclei", "n_germline_nuclei"):
        try:
            df[c] = pd.array(df[c].tolist(), dtype="Int64")      # whole counts, NA for a failed timepoint
        except (TypeError, ValueError):
            pass
    return df
