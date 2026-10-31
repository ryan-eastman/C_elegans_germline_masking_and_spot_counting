"""Batch hygiene (docs/ROADMAP_modular_pipeline.md step 5): file discovery with an exclusions file,
resume keyed on the completion marker, and the collect step that stacks per-image tables.

Shared by `germquant batch`, `germquant collect` and workflow/Snakefile so the three never disagree on
which images are in a study or on what "already done" means.
"""
from __future__ import annotations

import fnmatch
import json
import logging
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import provenance, schema
from .fsutil import long_path
from .stages import STAGES, stage_enabled

log = logging.getLogger(__name__)

# columns of batch_summary.csv (unchanged since the first release; R scripts read it)
SUMMARY_COLS = ["image_id", "n_nuclei", "n_germline", "qc_pass", "qc_flags", "out_dir"]
_SHORT_RE = re.compile(r"(noHS|HS)_(male|herm)_?(\d+)")
_BATCH_RE = re.compile(r"^(\d{8})_")


@dataclass
class Exclusions:
    """Study exclusions. `patterns` match a file stem only as WHOLE tokens: the whole stem, the parsed
    short id (``HS_male_07`` from ``..._HS_male_07``), or a contiguous run of ``_``-delimited tokens.
    So ``HS_herm_14`` never drops ``noHS_herm_14`` and ``HS_male_1`` never drops ``HS_male_10``.
    When `batch` is set (the analysis exclusions.json format), a stem is dropped only if its leading
    8-digit date equals `batch` as well (the same rule as analysis/coloc/chload.is_excluded)."""
    patterns: list[str] = field(default_factory=list)
    batch: str | None = None

    def matches(self, stem: str) -> str | None:
        """The pattern that excludes `stem`, or None."""
        if not self.patterns:
            return None
        if self.batch is not None:
            m = _BATCH_RE.match(stem)
            if not m or m.group(1) != str(self.batch):
                return None
        short = None
        m = _SHORT_RE.search(stem)
        if m:
            short = f"{m.group(1)}_{m.group(2)}_{m.group(3)}"
        tokens = stem.split("_")
        for pat in self.patterns:
            if pat == stem or (short is not None and pat == short):
                return pat
            pt = pat.split("_")
            n = len(pt)
            if n and any(tokens[i:i + n] == pt for i in range(len(tokens) - n + 1)):
                return pat
        return None


def load_exclusions(path: str | Path | None, base_dir: Path | None = None) -> Exclusions:
    """Read the study exclusions file: a JSON list of ids / token patterns, or an object with
    ``excluded_short_ids`` (and optionally ``excluded_batch``, the analysis/coloc/exclusions.json
    format). Missing path -> no exclusions."""
    if not path:
        return Exclusions()
    p = Path(path)
    if not p.is_absolute() and base_dir is not None:
        p = base_dir / p
    if not p.is_file():
        raise FileNotFoundError(f"qc.exclusions_file not found: {p}")
    data = json.loads(p.read_text(encoding="utf-8"))
    if isinstance(data, list):
        return Exclusions([str(x) for x in data])
    return Exclusions([str(x) for x in data.get("excluded_short_ids", [])],
                      str(data["excluded_batch"]) if data.get("excluded_batch") else None)


def discover_files(root: str | Path, glob: str = "**/*.nd2", exclude_patterns: list[str] | None = None,
                   exclusions: Exclusions | Iterable[str] | None = None, *, suffixes: Iterable[str] = (".nd2",),
                   skip_root: str | Path | None = None) -> tuple[list[Path], list[tuple[Path, str]]]:
    """Files to process, sorted, plus the (file, reason) pairs that were dropped. ``exclude_patterns``
    are filename globs (10x overviews, largeimage); ``exclusions`` are the study's Exclusions (a plain
    list of patterns is accepted too), matched against the image id (io.formats.image_stem).

    ``suffixes`` (io.input_suffixes, docs/DESIGN_system_packs.md 7.8) are the formats kept, matched on the whole
    format suffix in any letter case (io.formats.image_suffix: x.ome.tif is .ome.tif, never .tif). The default keeps
    exactly the .nd2 files of 8eb73de, under the same ids (image_stem is Path.stem for an .nd2). germquant's own
    output images (project.images.OUTPUT_SUFFIXES: <id>__nuclei_labels.tif and the like) are never kept; none of them
    ends in .nd2, so the default set is unchanged. With ``skip_root`` (the batch's output root, given only when
    io.input_suffixes is set, PLAN.md D20) nothing under that folder is kept either, when it lies inside ``root``
    (results in a folder of the input tree); an output root at or above ``root`` (results next to the images, or
    around them) holds the inputs themselves, so there only germquant's own images are left out.

    Two files of one folder whose image ids are equal (in any letter case) but whose formats differ (x.nd2 and an
    export x.tif) would share one results folder: the first in sorted order is kept, and each other one is dropped
    with the reason ``duplicate_id``. So is a later file of the same format whose id is equal in letter case too
    (x.tif and x.TIF, two files where the file system tells letter case apart), when the suffixes are not the
    default set. With the default suffixes every file is an .nd2 and neither rule applies, so a default batch
    keeps 8eb73de's files (PLAN.md D20)."""
    from .io.formats import image_stem, image_suffix
    from .project.images import is_output_name

    root = Path(root)
    excl = exclusions if isinstance(exclusions, Exclusions) else Exclusions([str(x) for x in (exclusions or [])])
    want = tuple(str(s).lower() for s in suffixes)
    skip = None if skip_root is None else os.path.normcase(os.path.abspath(skip_root))
    if skip is not None and not skip.startswith(os.path.normcase(os.path.abspath(root)).rstrip(os.sep) + os.sep):
        skip = None                        # at or above the input folder, or beside it: nothing of root lies under it
    keep: list[Path] = []
    dropped: list[tuple[Path, str]] = []
    first: dict[tuple[str, str], Path] = {}     # (folder, image id in any letter case) -> the file kept for it
    held: dict[tuple[str, str], Path] = {}      # (folder, image id) -> the file kept for it
    exact = want != (".nd2",)                   # io.input_suffixes set: one id twice in one format is a twin too
    for f in sorted(root.glob(glob)):
        if image_suffix(f) not in want or is_output_name(f.name):
            continue
        if skip is not None:
            full = os.path.normcase(os.path.abspath(f))
            if full == skip or full.startswith(skip.rstrip(os.sep) + os.sep):
                continue
        name = f.name.lower()
        pat = next((p for p in (exclude_patterns or []) if fnmatch.fnmatch(name, p.lower())), None)
        if pat is not None:
            dropped.append((f, f"exclude_pattern {pat}")); continue
        hit = excl.matches(image_stem(f))
        if hit is not None:
            dropped.append((f, f"exclusions_file {hit}")); continue
        here, iid = str(f.parent), image_stem(f)
        other = first.setdefault((here, iid.casefold()), f)
        twin = (other if other is not f and image_suffix(other) != image_suffix(f)
                else held.get((here, iid)) if exact else None)
        if twin is not None:
            # one image id in two formats (x.nd2 and x.tif), or twice in one format where the file system tells letter
            # case apart (x.tif and x.TIF): both would write the one results folder <id>/
            dropped.append((f, f"duplicate_id {iid}: {twin.name} in this folder has the same image id, "
                               f"which names one results folder")); continue
        held[(here, iid)] = f
        keep.append(f)
    return keep, dropped


def discover_options(cfg, out_root: str | Path | None = None) -> dict:
    """The keywords a recipe adds to discover_files, so `germquant batch` and workflow/Snakefile discover alike: {}
    when io.input_suffixes is unset (a default batch discovers the .nd2 files it always has, D20), else
    ``{"suffixes": (...), "skip_root": out_root}``. ValueError for a malformed io.input_suffixes
    (recipe_keys.input_suffixes_problem), before any file is discovered."""
    raw = cfg.get("io.input_suffixes")
    if raw is None:
        return {}
    from .recipe_keys import input_suffixes_problem

    problem = input_suffixes_problem(raw)
    if problem:
        raise ValueError("recipe refused before any file was discovered: " + problem)
    kw: dict = {"suffixes": tuple(raw)}
    if out_root is not None:
        kw["skip_root"] = Path(out_root)
    return kw


def enabled_stage_names(cfg) -> list[str]:
    return [s.name for s in STAGES if stage_enabled(cfg, s.name)]


def is_done(out_dir: str | Path, image_id: str, cfg, *, xy_stride: int = 1, z_range=None,
            segment_hash: str | None = None, retry_failed: bool = True, t: int | None = None,
            segment_hash_planar: str | None = None, sample_row_sha256: str | None = None) -> tuple[bool, str]:
    """True when ``<image_id>__done.json`` exists in out_dir with the same config_hash, the same set of
    enabled stages, the same run geometry (xy_stride, z_range), no failed stage (unless
    retry_failed=False) and, when given, the same segment sub-hash (model file or Cellpose version).
    ``t`` (a time-lapse timepoint, timelapse.series) additionally requires the marker's t_index to be it.

    A single-plane run's marker (extra ``"planar": true``, docs/PLANAR.md) carries the segment sub-hash of the
    planar payload, so it is compared with ``segment_hash_planar`` (None by default). A caller that gives only
    ``segment_hash`` (the time-lapse loop, an outside caller) has it compared with that hash, which a planar
    marker never matches: the image is reprocessed with the reason "planar flag", never skipped on a hash that
    was not checked. Without any segment hash the check is skipped, as before.

    ``sample_row_sha256`` (None by default; the image's sample sheet row, germquant.io.sample_sheet) is compared
    with the marker's: a sheet run records the row's sha256, or null for an image the sheet has no row for. The two
    must agree, None standing for no row, so a corrected design cell reruns that image, and so does a row added for
    an image that ran without one, or a row (or the whole sheet) taken away from an image that ran with one; the
    reason starts with "sample row". A marker of a run without a sheet, checked without a row, is unchanged."""
    m = provenance.read_done_marker(out_dir, image_id)
    if not m:
        return False, "no completion marker"
    if m.get("config_hash") != cfg.hash:
        return False, f"config_hash {m.get('config_hash')} != {cfg.hash}"
    if sorted(m.get("enabled_stages", [])) != sorted(enabled_stage_names(cfg)):
        return False, "enabled stages differ"
    if "xy_stride" not in m or int(m["xy_stride"]) != int(xy_stride):
        return False, f"xy_stride {m.get('xy_stride', 'unknown')} != {xy_stride}"
    want = list(z_range) if z_range else None
    if (m.get("z_range") or None) != want:
        return False, "z_range differs"
    if t is not None and m.get("t_index") != int(t):
        return False, f"t_index {m.get('t_index', 'absent')} != {int(t)}"
    if retry_failed and m.get("failed_stages"):
        return False, f"stage(s) failed last run: {', '.join(m['failed_stages'])}"
    if sample_row_sha256 is not None or m.get("sample_row_sha256") is not None:
        # the sample sheet row (docs/PROJECTS.md, "Sample sheet"): a marker without the key ran without a sheet
        why = _sample_row_differs(m, sample_row_sha256)
        if why:
            return False, why
    if m.get("planar") is True and (segment_hash is not None or segment_hash_planar is not None):
        # a single-plane marker: its segment sub-hash is the planar one (the safe default of D8 without it)
        want = segment_hash_planar if segment_hash_planar is not None else segment_hash
        if m.get("stage_hashes", {}).get("segment") != want:
            if segment_hash_planar is None:
                return False, "planar flag: a single-plane marker, and no planar segment hash was given to compare"
            return False, _segment_differs(cfg)
        return True, "done"
    if segment_hash is not None and m.get("stage_hashes", {}).get("segment") != segment_hash:
        return False, _segment_differs(cfg)
    return True, "done"


def _segment_differs(cfg) -> str:
    """The reason is_done gives for a segment sub-hash that differs. A labels run (segmentation.nuclei.labels_from,
    docs/LABELS.md) resolves no model, so there the labels file or the segment stage's tools moved it (D52); every
    other run keeps the text it always had."""
    if cfg.get("segmentation.nuclei.labels_from") is not None:
        return "segment stage hash differs (the labels file, or the Cellpose or torch version, changed)"
    return "segment stage hash differs (model file or Cellpose version changed)"


def _sample_row_differs(marker: dict, given: str | None) -> str | None:
    """Why the marker's sample sheet row is not the ``given`` one (None: the image has no row now), or None."""
    have = marker.get("sample_row_sha256")
    if have == given:
        return None
    if given is None:
        return ("sample row: the image ran with a sample sheet row, and none was given to compare (no sample sheet, "
                "or no row for this image in it)")
    if "sample_row_sha256" not in marker:
        return "sample row differs (the image ran without a sample sheet)"
    if have is None:
        return "sample row differs (the image had no row in the sample sheet when it ran)"
    return "sample row differs (this image's row in the sample sheet changed)"


def apply_framing_qc(rows: list[dict], factor: float = 1.8) -> list[dict]:
    """Batch-level framing QC: flag images whose germline count is a strong outlier vs the batch median
    (two gonad arms, extra tissue, fuller distal capture: worth eyeballing the montage; spot counts are
    unaffected, so advisory). The same rule for `batch`, `batch --resume` and `collect`."""
    germ = [r["n_germline"] for r in rows if r.get("qc_pass") and r.get("n_germline", 0) > 0]
    if len(germ) >= 4:
        import statistics

        med = statistics.median(germ)
        for r in rows:
            if med > 0 and r.get("n_germline", 0) > factor * med:
                flag = f"qc:germline_count_outlier_{r['n_germline']}_vs_median{med:.0f}_review_framing"
                r["qc_flags"] = f"{r['qc_flags']};{flag}" if r.get("qc_flags") else flag
    return rows


def _find(out_root: Path, suffix: str) -> list[Path]:
    """Files ending in `suffix` under out_root, walked through the extended-length path so deep NAS
    trees beyond MAX_PATH do not fail; results are ordinary paths."""
    root_lp = long_path(out_root)
    hits: list[Path] = []
    for dirpath, _dirs, names in os.walk(root_lp):
        rel = os.path.relpath(dirpath, root_lp)
        base = out_root if rel == "." else out_root / rel
        for n in names:
            if n.endswith(suffix) and not n.startswith("batch_"):
                hits.append(base / n)
    return sorted(hits)


def _int(v, default=0) -> int:
    try:
        return int(v) if pd.notna(v) else default
    except (TypeError, ValueError):
        return default


def _bool(v) -> bool:
    if isinstance(v, str):
        return v.strip().lower() == "true"
    try:
        return bool(v) if pd.notna(v) else False
    except (TypeError, ValueError):
        return False


def collect_series_summaries(out_root: str | Path) -> int | None:
    """Stack every ``<stem>__series_summary.csv`` under out_root (one per time-lapse series,
    timelapse.series.process_series) into ``batch_series_summary.csv``. Returns the row count, or None
    when the tree holds no series (then nothing is written, so a fixed-image study's folder is unchanged)."""
    out_root = Path(out_root)
    files = _find(out_root, "__series_summary.csv")
    if not files:
        return None
    frames = []
    for f in files:
        try:
            frames.append(pd.read_csv(long_path(f), low_memory=False))
        except Exception as e:  # noqa: BLE001 - one unreadable table must not block the stack
            log.warning("collect: skipping %s (%s)", f, e)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=schema.SERIES_SUMMARY)
    df.to_csv(long_path(out_root / "batch_series_summary.csv"), index=False)
    return len(df)


def _join_sheet(df: pd.DataFrame, table: str, frames: list, frame_ids: list[str], columns: list[str],
                values: dict[str, dict[str, str]]) -> pd.DataFrame:
    """batch_<table>.csv with a sample sheet's extra ``columns`` appended: every row takes the value its image's row
    gives (``values``, by the image id of the table file the row came from), "" for an image without a row. A
    column the table already has is left as it is, with a warning."""
    ids = [iid for iid, frame in zip(frame_ids, frames) for _ in range(len(frame))]
    new = {}
    for c in columns:
        if c in df.columns:
            log.warning("collect: the sample sheet's column %s is already a column of batch_%s.csv; it is not joined "
                        "there", c, table)
            continue
        new[c] = [values.get(i, {}).get(c, "") for i in ids]
    if not new:
        return df
    return pd.concat([df, pd.DataFrame(new, index=df.index, dtype=object)], axis=1)


def collect(out_root: str | Path, tables: tuple[str, ...] = tuple(schema.TABLES), *,
            image_ids: set[str] | None = None, extra_rows: list[dict] | None = None,
            outlier_factor: float | None = 1.8, require_marker: bool = True, sample_sheet=None) -> dict[str, int]:
    """Stack every per-image ``<image_id>__<table>.csv`` under out_root into ``batch_<table>.csv`` and
    rebuild ``batch_summary.csv`` (image_id, n_nuclei, n_germline, qc_pass, qc_flags, out_dir).

    * only folders holding a completion marker are stacked (``require_marker``); a folder with tables
      but no marker gets a summary row with qc_pass False and flag ``incomplete:no_done_marker``;
    * ``image_ids`` restricts the collection to a study (the files the batch discovered, so images
      excluded after a first pass drop out);
    * ``extra_rows`` (this pass's EXCEPTION rows) override or add summary rows by image_id;
    * the framing outlier QC is applied over the final rows (``outlier_factor``; None = skip);
    * time-lapse series (timelapse.series): the per-timepoint folders stack like any image (their ids
      carry the _t<NNN> suffix), and when any ``<stem>__series_summary.csv`` is present those are stacked
      into ``batch_series_summary.csv`` as well (counted under ``series_summary``; nothing is written
      without one);
    * ``sample_sheet`` (None by default; a SampleSheet or a sheet's path, germquant.io.sample_sheet): its extra
      columns (neither the key nor a design column) are appended to every ``batch_<table>.csv``, joined by the
      image id of each table file, empty for an image without a row. A column a stacked table already has is not
      joined there (a warning says so); the design columns stay the ones each image ran with, and
      batch_summary.csv keeps its six columns.
    Idempotent; safe to rerun after any subset of images was reprocessed. Returns {table: n_rows}."""
    out_root = Path(out_root)
    sheet_columns, sheet_values = [], {}
    if sample_sheet is not None:
        from .io.sample_sheet import SampleSheet, read_sample_sheet

        if not isinstance(sample_sheet, SampleSheet):
            sample_sheet = read_sample_sheet(sample_sheet)
        sheet_columns, sheet_values = sample_sheet.extra_columns, sample_sheet.extra_values()
    counts: dict[str, int] = {}
    summaries = _find(out_root, "__image_summary.csv")
    folders: dict[str, tuple[Path, dict | None]] = {}
    for f in summaries:
        iid = f.name.split("__")[0]
        if image_ids is not None and iid not in image_ids:
            continue
        folders[iid] = (f.parent, provenance.read_done_marker(f.parent, iid))
    complete = {iid for iid, (_d, m) in folders.items() if m or not require_marker}
    # extra spot-instance tables (spots_<name>) are discovered from disk so no schema change is needed
    extra = sorted({f.name.split("__", 1)[1][:-4] for f in _find(out_root, ".csv")
                    if "__spots_" in f.name and f.name.split("__")[0] in complete})
    for t in tuple(tables) + tuple(x for x in extra if x not in tables):
        frames = []
        frame_ids: list[str] = []          # the image id of each frame (its file name), for the sample sheet join
        for f in _find(out_root, f"__{t}.csv"):
            iid = f.name.split("__")[0]
            if iid not in complete:
                continue
            try:
                frames.append(pd.read_csv(long_path(f), low_memory=False))
                frame_ids.append(iid)
            except Exception as e:  # noqa: BLE001 - one unreadable table must not block the stack
                log.warning("collect: skipping %s (%s)", f, e)
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=schema.TABLES.get(t, schema.SPOTS))
        if sheet_columns:
            df = _join_sheet(df, t, frames, frame_ids, sheet_columns, sheet_values)
        df.to_csv(long_path(out_root / f"batch_{t}.csv"), index=False)
        counts[t] = len(df)
    n_series = collect_series_summaries(out_root)
    if n_series is not None:
        counts["series_summary"] = n_series
    rows: dict[str, dict] = {}
    for iid, (d, marker) in folders.items():
        f = d / f"{iid}__image_summary.csv"
        try:
            s = pd.read_csv(long_path(f), low_memory=False)
            r = s.iloc[0] if not s.empty else {}
        except Exception as e:  # noqa: BLE001 - keep stacking the rest
            log.warning("collect: skipping %s (%s)", f, e)
            continue
        flags = r.get("qc_flags", "") if pd.notna(r.get("qc_flags", "")) else ""
        row = {"image_id": r.get("image_id", iid) if pd.notna(r.get("image_id", iid)) else iid,
               "n_nuclei": _int(r.get("n_nuclei", 0)), "n_germline": _int(r.get("n_germline_nuclei", 0)),
               "qc_pass": _bool(r.get("qc_pass", False)), "qc_flags": str(flags), "out_dir": str(d)}
        if not marker and require_marker:
            row["qc_pass"] = False
            row["qc_flags"] = f"{row['qc_flags']};incomplete:no_done_marker" if row["qc_flags"] else "incomplete:no_done_marker"
        rows[str(row["image_id"])] = row
    for r in extra_rows or []:
        rows[str(r["image_id"])] = {k: r.get(k, "") for k in SUMMARY_COLS}
    out_rows = list(rows.values())
    if outlier_factor is not None:
        out_rows = apply_framing_qc(out_rows, outlier_factor)
    pd.DataFrame(out_rows, columns=SUMMARY_COLS).to_csv(long_path(out_root / "batch_summary.csv"), index=False)
    counts["batch_summary"] = len(out_rows)
    return counts
