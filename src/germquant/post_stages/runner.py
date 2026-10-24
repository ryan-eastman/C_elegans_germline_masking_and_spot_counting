"""Running post-stages on a finished result folder (docs/POST_STAGES.md; design 5.6; decisions D28 and D32 to D34).

``run_post_stages(out_dir, image_id, entries, *, cfg, image_path, stack=None)`` runs each listed post-stage on what
process_image wrote into ``out_dir``: the tables and the label image, and the image itself read through the run's
reader at the run's xy_stride and z_range (``stack`` hands over the stack a caller already holds, through
process_image's on_stack hook, D34). A post-stage is ``fn(inputs, params) -> {table: DataFrame}``; the engine,
not the post-stage, writes each returned table as <image_id>__post_<name>_<table>.csv with exactly the declared
columns, and records every stage in <image_id>__post.json:

    {"format": 1, "image_id", "post_stages_file_sha256",
     "stages": {name: {"version", "params", "sub_hash", "done", "error", "tables", "finished", "inputs"}}}

The sub-hash (12 hex, like the stage sub-hashes) covers the post-stage's version and declaration, its resolved
parameters, the post_stages.yaml digest and the digests of its inputs (D32): the main run's read sub-hash, the
spacing and the channel roles it reads, the LF-normalised sha256 of each table it declares and of the label image
when it declares it, and the run's stride and z range when it reads the image. The raw image is not hashed
again. Under resume (the default) a stage whose recorded sub-hash matches, that is done and whose tables are on
disk, is skipped; anything else runs again.

What this never does: edit a file process_image wrote (the tables, the label images, the done marker, the stage
record), fail the image, or raise for a post-stage's error. A failure (a post-stage that is not installed, a
parameter out of range, an exception, a table outside the declaration, a post-stage that writes into the folder
itself) is recorded in __post.json with done false; the other stages still run. The runner's own guards refuse a
folder without a completion marker, a time-series timepoint (refused in T5, D33) and a folder made with another
config_hash, recording the refusal for every stage.

Imports: run_post_stages imports germquant.pipeline before anything else (D28), so torch claims its cuBLAS before
a reader, CuPy or the post-stage's own module can. This module itself imports nothing heavy.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from pathlib import Path

from .. import provenance
from ..fsutil import long_path
from .config import (
    as_declarations,
    canonical_text,
    entries_problems,
    file_sha256,
    load_post_stages_file,
    problems_by_name,
    text_sha256,
    write_text_atomic,
)
from .declaration import resolve_params
from .discovery import PostStageError, load_post_stage

log = logging.getLogger(__name__)

RECORD_FORMAT = 1
LABELS_FILE = "{image_id}__nuclei_labels.tif"
SERIES_REFUSAL = ("post-stages are refused for time series in this release (design 5.6): run the series without "
                  "post-stages")


class PostStageContractError(RuntimeError):
    """A post-stage that reads what it did not declare, returns what it did not declare, or changes the results
    folder itself."""


# ---- names and records -------------------------------------------------------------------------------------

def post_record_path(out_dir, image_id: str) -> Path:
    return Path(out_dir) / f"{image_id}__post.json"


def post_table_name(image_id: str, name: str, table: str) -> str:
    """The file name of a post table: <image_id>__post_<name>_<table>.csv."""
    return f"{image_id}__post_{name}_{table}.csv"


def read_post_record(out_dir, image_id: str) -> dict | None:
    """The parsed <image_id>__post.json, or None when there is none or it cannot be read."""
    try:
        with open(long_path(post_record_path(out_dir, image_id)), encoding="utf-8") as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        return None
    return rec if isinstance(rec, dict) else None


def sub_hash(payload: dict) -> str:
    """sha256 (12 hex) of the canonical JSON of a sub-hash payload."""
    text = json.dumps(payload, sort_keys=True, ensure_ascii=True, separators=(",", ":"), default=str)
    return hashlib.sha256(text.encode("ascii")).hexdigest()[:12]


# ---- inputs ------------------------------------------------------------------------------------------------

def _table_path(out_dir: Path, image_id: str, table: str) -> Path:
    csv = out_dir / f"{image_id}__{table}.csv"
    if not os.path.isfile(long_path(csv)):
        parquet = out_dir / f"{image_id}__{table}.parquet"
        if os.path.isfile(long_path(parquet)):
            return parquet
    return csv


def _digest(path: Path) -> str:
    """sha256 of a file, LF-normalised for a CSV table; "absent" when there is no file."""
    try:
        data = Path(long_path(path)).read_bytes()
    except FileNotFoundError:
        return "absent"
    if path.suffix.lower() == ".csv":
        data = data.replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def _read_table(path: Path):
    import pandas as pd

    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(long_path(path))
    return pd.read_csv(long_path(path), low_memory=False, float_precision="round_trip")


def _run_facts(out_dir: Path, image_id: str) -> tuple[tuple[float, float, float] | None, dict]:
    """(spacing (dz, dy, dx) in um, {role: channel index or None}) as the main run recorded them in its
    image_summary table (the voxel_d*_um and channel_map_json columns), so every route reads the same values."""
    path = _table_path(out_dir, image_id, "image_summary")
    try:
        df = _read_table(path)
    except Exception:  # noqa: BLE001 - a missing or damaged summary leaves the facts unknown
        return None, {}
    if df.empty:
        return None, {}
    row = df.iloc[0]
    spacing = None
    try:
        vals = [float(row[c]) for c in ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um")]
        if all(v == v for v in vals):              # not NaN
            spacing = (vals[0], vals[1], vals[2])
    except (KeyError, TypeError, ValueError):
        spacing = None
    roles: dict = {}
    try:
        parsed = json.loads(str(row["channel_map_json"]))
        if isinstance(parsed, dict):
            roles = {str(k): (None if v is None else int(v)) for k, v in parsed.items()}
    except (KeyError, TypeError, ValueError):
        roles = {}
    return spacing, roles


def _read_only(arr):
    v = arr.view()
    v.flags.writeable = False
    return v


class PostInputs:
    """What a post-stage receives (docs/POST_STAGES.md):

    * ``image_id``, ``out_dir`` (a Path), ``spacing`` ((dz, dy, dx) in um as the run recorded it, or None) and
      ``roles`` ({role: channel index or None}, the run's resolved channel map);
    * ``tables(name)``: a declared table of the main run as a DataFrame (a fresh copy each call);
    * ``labels()``: the nucleus label image as a read-only int32 (Z, Y, X) array (a single plane as (1, Y, X));
    * ``image(role)``: the channel of a declared role as a read-only (Z, Y, X) array, read through the run's
      reader at the run's xy_stride and z_range (once per image, however many post-stages read it).

    Reading anything the declaration does not name raises PostStageContractError."""

    def __init__(self, *, image_id: str, out_dir: Path, spacing, roles: dict, reads: dict, cache: dict,
                 image_loader):
        self.image_id = image_id
        self.out_dir = Path(out_dir)
        self.spacing = tuple(spacing) if spacing is not None else None
        self.roles = dict(roles)
        self._reads = reads
        self._cache = cache
        self._image_loader = image_loader

    def tables(self, name: str):
        if name not in self._reads.get("tables", []):
            raise PostStageContractError(f"table {name!r} is not declared in reads.tables "
                                         f"({self._reads.get('tables', [])})")
        key = ("table", name)
        if key not in self._cache:
            path = _table_path(self.out_dir, self.image_id, name)
            if not os.path.isfile(long_path(path)):
                raise FileNotFoundError(f"the main run wrote no {name} table ({path.name})")
            self._cache[key] = _read_table(path)
        return self._cache[key].copy()

    def labels(self):
        if not self._reads.get("labels"):
            raise PostStageContractError("the label image is not declared (reads.labels is not true)")
        key = ("labels",)
        if key not in self._cache:
            import numpy as np
            import tifffile

            path = self.out_dir / LABELS_FILE.format(image_id=self.image_id)
            if not os.path.isfile(long_path(path)):
                raise FileNotFoundError(f"the main run wrote no label image ({path.name})")
            arr = np.asarray(tifffile.imread(long_path(path)))
            if arr.ndim == 2:
                arr = arr[None]
            self._cache[key] = arr.astype(np.int32, copy=False)
        return _read_only(self._cache[key])

    def image(self, role: str):
        if not self._reads.get("image") or role not in self._reads.get("roles", []):
            raise PostStageContractError(f"channel role {role!r} is not declared in reads.roles "
                                         f"({self._reads.get('roles', [])}) with reads.image true")
        return _read_only(self._image_loader(role))


def _image_loader(*, image_path, marker: dict, stack, roles: dict, cfg, cache: dict):
    def load(role: str):
        import numpy as np

        if ("image",) not in cache:
            if stack is not None:
                cache[("image",)] = (stack, list(getattr(stack, "channel_names", []) or []))
            else:
                if image_path is None:
                    raise FileNotFoundError("no image path was given, so the image cannot be read")
                from .. import pipeline

                zr = marker.get("z_range")
                zr = (int(zr[0]), int(zr[1])) if zr else None
                meta, st = pipeline._read_input(Path(image_path), xy_stride=int(marker.get("xy_stride", 1) or 1),
                                                z_range=zr)
                cache[("image",)] = (st, list(meta.get("channel_names") or getattr(st, "channel_names", [])))
        st, names = cache[("image",)]
        idx = roles.get(role)
        if role not in roles and cfg is not None and names:
            idx = cfg.channel_map.resolve(names)[0].get(role)
        if idx is None:
            raise PostStageContractError(f"channel role {role!r} is not bound in this run (roles {roles})")
        return np.asarray(st.data[int(idx)])

    return load


def input_digests(out_dir: Path, image_id: str, declaration: dict, marker: dict, spacing, roles: dict) -> dict:
    """What the sub-hash takes from the inputs (D32): the main run's read sub-hash, the spacing, the declared
    roles' channel indices, the digest of each declared table and of the label image when declared, and the
    stride and z range when the image is read. The raw image is not hashed."""
    reads = declaration["reads"]
    out = {"read_sub_hash": (marker.get("stage_hashes") or {}).get("read"),
           "spacing": list(spacing) if spacing is not None else None,
           "tables": {t: _digest(_table_path(out_dir, image_id, t)) for t in reads["tables"]}}
    if reads["roles"]:
        out["roles"] = {r: roles.get(r) for r in reads["roles"]}
    if reads["labels"]:
        out["labels"] = _digest(out_dir / LABELS_FILE.format(image_id=image_id))
    if reads["image"]:
        out["image"] = {"xy_stride": int(marker.get("xy_stride", 1) or 1),
                        "z_range": list(marker["z_range"]) if marker.get("z_range") else None}
    return out


# ---- one stage ---------------------------------------------------------------------------------------------

def _listing(folder: Path) -> dict[str, tuple[int, int]]:
    out = {}
    with os.scandir(long_path(folder)) as it:
        for e in it:
            if e.is_file():
                st = e.stat()
                out[e.name] = (st.st_size, st.st_mtime_ns)
    return out


def _conform(name: str, declaration: dict, produced) -> dict:
    """The returned tables, checked against the declaration and put in the declared column order."""
    import pandas as pd

    if not isinstance(produced, dict):
        raise PostStageContractError(f"post-stage {name} must return a mapping {{table: DataFrame}}, got "
                                     f"{type(produced).__name__}")
    declared = declaration["writes"]
    missing = [t for t in declared if t not in produced]
    if missing:
        raise PostStageContractError(f"post-stage {name} did not return the declared table(s) {missing}")
    extra = sorted(str(t) for t in produced if t not in declared)
    if extra:
        raise PostStageContractError(f"post-stage {name} returned table(s) {extra} that writes does not declare")
    out = {}
    for table, cols in declared.items():
        df = produced[table]
        if not isinstance(df, pd.DataFrame):
            raise PostStageContractError(f"post-stage {name}: table {table} must be a pandas DataFrame, got "
                                         f"{type(df).__name__}")
        have = [str(c) for c in df.columns]
        if len(set(have)) != len(have):
            raise PostStageContractError(f"post-stage {name}: table {table} holds a column twice")
        lacking = [c for c in cols if c not in have]
        if lacking:
            raise PostStageContractError(f"post-stage {name}: table {table} lacks the declared column(s) {lacking}")
        undeclared = [c for c in have if c not in cols]
        if undeclared:
            raise PostStageContractError(f"post-stage {name}: table {table} holds column(s) {undeclared} that "
                                         f"writes.{table} does not declare")
        out[table] = df[list(cols)]
    return out


def _write_tables(out_dir: Path, image_id: str, name: str, frames: dict) -> list[str]:
    """Write every table to a temporary file first, then move each into place, so a table that cannot be
    written leaves no post table of this stage changed."""
    staged = []
    try:
        for table, df in frames.items():
            final = out_dir / post_table_name(image_id, name, table)
            tmp = final.with_name(final.name + ".part")
            df.to_csv(long_path(tmp), index=False)
            staged.append((tmp, final))
        for tmp, final in staged:
            os.replace(long_path(tmp), long_path(final))
    finally:
        for tmp, _final in staged:
            try:
                os.remove(long_path(tmp))
            except FileNotFoundError:
                pass
    return [final.name for _tmp, final in staged]


def _run_one(name: str, fn, declaration: dict, params: dict, inputs: PostInputs, out_dir: Path,
             image_id: str) -> tuple[bool, str | None, list[str]]:
    """(done, error, tables written) of one post-stage call."""
    before = _listing(out_dir)
    try:
        produced = fn(inputs, dict(params))
        after = _listing(out_dir)
        if after != before:
            changed = sorted(set(after) ^ set(before) | {k for k in before if k in after and after[k] != before[k]})
            raise PostStageContractError(f"post-stage {name} changed the results folder itself ({changed}): a "
                                         f"post-stage returns its tables and the engine writes them")
        frames = _conform(name, declaration, produced)
        return True, None, _write_tables(out_dir, image_id, name, frames)
    except Exception as e:  # noqa: BLE001 - a post-stage failure is recorded, never raised (design 5.6)
        log.warning("%s: post-stage %s failed (%s: %s)", image_id, name, type(e).__name__, e)
        return False, f"{type(e).__name__}: {e}", []


# ---- the runner --------------------------------------------------------------------------------------------

def _main_problem(marker: dict | None, cfg) -> str | None:
    if marker is None:
        return "no completion marker: the main run of this image has not finished, so no post-stage runs on it"
    if "t_index" in marker:
        return f"{SERIES_REFUSAL} (this folder holds timepoint {marker.get('t_index')})"
    want = getattr(cfg, "hash", None) if cfg is not None else None
    if want is not None and marker.get("config_hash") != want:
        return (f"the result folder was made with config_hash {marker.get('config_hash')}, not with this run's "
                f"config_hash {want}")
    return None


def run_post_stages(out_dir, image_id: str, entries, *, cfg, image_path, stack=None, resume: bool = True) -> dict:
    """Run the post-stages ``entries`` on the finished result folder ``out_dir`` of image ``image_id``.

    ``entries`` is a post-stage list ([{"name", "params"}]) or the path of a post_stages.yaml (whose file digest
    then enters the sub-hashes; a list digests like the file write_post_stages_file makes of it). ``cfg`` is the
    run's Config (its config_hash must be the marker's), ``image_path`` the image (read only when a post-stage
    reads it and no ``stack`` is handed over). ``resume`` (default True) skips a stage whose record is done
    with the same sub-hash and whose tables exist.

    Returns {"image_id", "record" (the path of __post.json), "ran", "skipped", "failed" (post-stage names),
    "post" (the record, None when nothing was listed)}, plus "error" when the list itself cannot be read. An
    empty list writes nothing and calls nothing. Never raises for a post-stage's failure."""
    from .. import pipeline  # noqa: F401  (D28: torch first, before any reader, CuPy or post-stage module)

    out_dir = Path(out_dir)
    result = {"image_id": image_id, "record": str(post_record_path(out_dir, image_id)), "ran": [], "skipped": [],
              "failed": [], "post": None}
    file_sha = None
    if entries is None:
        entries = []
    elif isinstance(entries, (str, os.PathLike)):
        try:
            source = Path(entries)
            entries = load_post_stages_file(source)
            file_sha = file_sha256(source)
        except (OSError, ValueError) as e:
            log.warning("%s: post-stages not run: %s", image_id, e)
            result["error"] = f"{type(e).__name__}: {e}"
            return result
    entries = list(entries)
    if not entries:
        return result
    shape = entries_problems(entries)
    if shape:
        log.warning("%s: post-stages not run: %s", image_id, "; ".join(shape))
        result["error"] = "; ".join(shape)
        return result
    decls, reasons = as_declarations(None)
    if file_sha is None:
        file_sha = text_sha256(canonical_text([{"name": e["name"], "params": dict(e.get("params") or {})}
                                               for e in entries], decls))
    problems: dict[str, list[str]] = {}
    for name, prob in problems_by_name(entries, decls, reasons):
        problems.setdefault(name, []).append(prob)
    old = read_post_record(out_dir, image_id) or {}
    old_stages = old.get("stages") if old.get("format") == RECORD_FORMAT and isinstance(old.get("stages"), dict) else {}
    marker = provenance.read_done_marker(out_dir, image_id)
    refused = _main_problem(marker, cfg)
    spacing, roles = _run_facts(out_dir, image_id)
    if spacing is None and stack is not None and getattr(stack, "spacing", None) is not None:
        spacing = tuple(float(s) for s in stack.spacing)
    cache: dict = {}
    loader = _image_loader(image_path=image_path, marker=marker or {}, stack=stack, roles=roles, cfg=cfg,
                           cache=cache)
    stages: dict[str, dict] = {}
    for entry in entries:
        name = entry["name"]
        if name in stages:
            continue                                   # listed twice: recorded once, as a failure
        given = dict(entry.get("params") or {})
        decl = decls.get(name)
        rec = {"version": decl["version"] if decl else None, "params": given, "sub_hash": None, "done": False,
               "error": None, "tables": [], "finished": None}
        why = refused or ("; ".join(problems[name]) if name in problems else None)
        if why is None:
            params = resolve_params(decl, given)
            inputs = input_digests(out_dir, image_id, decl, marker, spacing, roles)
            sh = sub_hash({"name": name, "version": decl["version"], "reads": decl["reads"],
                           "writes": decl["writes"], "params": params, "post_stages_file_sha256": file_sha,
                           "inputs": inputs})
            rec.update(params=params, sub_hash=sh, inputs=inputs)
            prev = old_stages.get(name)
            if (resume and isinstance(prev, dict) and prev.get("done") is True and prev.get("sub_hash") == sh
                    and prev.get("tables")
                    and all(os.path.isfile(long_path(out_dir / t)) for t in prev["tables"])):
                stages[name] = prev
                result["skipped"].append(name)
                log.info("%s: post-stage %s is up to date (sub-hash %s)", image_id, name, sh)
                continue
            t0 = time.perf_counter()
            try:
                fn, decl = load_post_stage(name)
            except PostStageError as e:
                done, error, tables = False, f"{type(e).__name__}: {e}", []
            else:
                post_inputs = PostInputs(image_id=image_id, out_dir=out_dir, spacing=spacing, roles=roles,
                                         reads=decl["reads"], cache=cache, image_loader=loader)
                done, error, tables = _run_one(name, fn, decl, params, post_inputs, out_dir, image_id)
            rec.update(done=done, error=error, tables=tables)
            log.info("%s: post-stage %s %s in %.2f s", image_id, name, "done" if done else "failed",
                     time.perf_counter() - t0)
        else:
            rec["error"] = why
            log.warning("%s: post-stage %s not run: %s", image_id, name, why)
        rec["finished"] = provenance.run_timestamp()
        stages[name] = rec
        result["ran" if rec["done"] else "failed"].append(name)
    record = {"format": RECORD_FORMAT, "image_id": image_id, "post_stages_file_sha256": file_sha, "stages": stages}
    if record != old:
        try:
            write_text_atomic(post_record_path(out_dir, image_id),
                              json.dumps(record, indent=1, ensure_ascii=True, default=str) + "\n")
        except Exception as e:  # noqa: BLE001 - the record is the post-stages' own; the image stays done
            log.warning("%s: could not write the post-stage record (%s: %s)", image_id, type(e).__name__, e)
            result["error"] = f"the post-stage record could not be written ({type(e).__name__}: {e})"
    result["post"] = record
    return result
