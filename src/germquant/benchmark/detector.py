"""Detector bench harness: run a spot detection backend on a cached gonad exactly as the pipeline's
spots stage ran it, and score the result against the Imaris truth.

Why a cache. The spots stage sees three things: the foci channel of the stack, the nucleus label image
and the voxel size. Re-running a whole gonad through the pipeline to try one detector setting costs a
Cellpose segmentation; the cache keeps exactly those inputs, taken from the frozen golden run, so a
detector can be re-run on them in minutes and its output is directly comparable with the deployed one.

The cache (``scripts/bench_cache.py`` writes it, ``load_cache`` reads it), one folder per gonad under
the cache root (``C:/Users/ryane/germquant_golden/bench_cache/<tag>/`` on the workstation):

* ``foci.npy`` and ``dna.npy``: the foci and DNA channels, (Z, Y, X) at full resolution in the dtype the
  reader returned (uint16 for the lab's .nd2 stacks), read with germquant's own reader. The channel
  indices come from the golden run's ``channel_map_json`` (its image summary), never from a constant.
* ``labels.npy``: int32 (Z, Y, X), exactly the label image the golden run handed to ``detect_spots``.
  In the pipeline (src/germquant/pipeline.py, the spots stage, at the pinned commit 8666f07 and today)
  that is the full nucleus segmentation, every segmented nucleus, germline or not: the germline stage
  only flags rows of the nuclei table (``in_germline``) and never edits the label image, and the same
  array is what the run saved as ``<image>__nuclei_labels.tif``. So the cache holds that TIFF as it is,
  and the germline ids are recorded next to it for scoring.
* ``nuclei.csv``: nucleus_id, in_germline, the centroid columns and the golden ``n_spots`` of every
  nucleus of the golden nuclei table (small; the scoring needs the germline flag and the centroids).
* ``meta.json``: image id, stem, voxel size, channel role indices, the marker name, the golden run's
  ``spots:`` config section (from its run_manifest.json), the .nd2 path, size and sha256 (from the
  pre-registered split record), the germline label ids, the shape, dtype and a sha256 of the data of
  every array, and where everything came from. It is written last, so a folder without it is an
  unfinished cache and ``load_cache`` refuses it.

Running (``run_backend``): ``germquant.spots.detect_spots`` with the cached foci channel, labels and
voxel size, the marker, the keyword parameters built from the cached ``spots:`` section with the same
expressions and code defaults as the pipeline's spots stage, ``backend`` and ``cfg_section``. ``params``
override keys of the section first (a nested mapping such as ``log:`` merges one level deep, the rule of
``pipeline.spot_backend_section``), so ``params={}`` is the golden configuration.

Scoring (``score``): ``germquant.benchmark.spots.benchmark_gonad`` with the detector's spots and counts in
place of the golden ones, against the truth files the manifest lists: Hungarian matching inside the
0.3 um (xy) by 0.4 um (z) cylinder, per-nucleus agreement (expert spots mapped into our labels against
our counts) over the germline nuclei and over the traced subset, with the truth counts taken from Class
A+B spots and from Class A spots only. ``pool_scores`` adds the cluster bootstrap by gonad
(``germquant.benchmark.stats``).

Test guard. Under docs/BENCHMARK_PROTOCOL.md the test gonads are scored once, by a frozen method. Every
entry point here refuses a test gonad unless the caller passes ``allow_test=True``, which only
``scripts/bench_detector.py --final`` does. The pre-registered roles come from
``germquant.validate.benchmark_split.PREREGISTERED_ROLES`` (the split never changes) and the manifest;
a disagreement between the two is an error.
"""
from __future__ import annotations

import copy
import hashlib
import json
import logging
import math
import os
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import spots as bs
from .manifest import Entry, PipelineOutputs, load_labels, load_manifest, load_truth
from .stats import agreement_metrics

log = logging.getLogger(__name__)

DEFAULT_CACHE_ROOT = Path("C:/Users/ryane/germquant_golden/bench_cache")
DEFAULT_MANIFEST = Path("C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/manifest.json")
CACHE_FORMAT = 1
CACHE_FILES = {"foci": "foci.npy", "dna": "dna.npy", "labels": "labels.npy", "nuclei": "nuclei.csv",
               "meta": "meta.json"}
ARRAYS = ("foci", "dna", "labels")
DEFAULT_Z_CHUNK = 8
# the pipeline output files the cache is built from (the names of germquant.validate.benchmark_split)
GOLDEN_FILES = {
    "nuclei": "{stem}__nuclei.csv",
    "nuclei_labels": "{stem}__nuclei_labels.tif",
    "spots": "{stem}__spots.csv",
    "spots_parquet": "{stem}__spots.parquet",
    "image_summary": "{stem}__image_summary.csv",
    "run_manifest": "run_manifest.json",
}
REQUIRED_GOLDEN = ("nuclei", "nuclei_labels", "spots", "image_summary", "run_manifest")
NUCLEI_COLUMNS = ["nucleus_id", "in_germline", "centroid_x_um", "centroid_y_um", "centroid_z_um", "n_spots"]
# docs/BENCHMARK_PROTOCOL.md section 6.3: 1000 resamples by gonad, seed 20251105
N_BOOT = 1000
BOOT_SEED = 20251105
# the per-spot columns compared with the golden table by compare_with_golden
COMPARE_COLUMNS = ["spot_id", "nucleus_id", "marker", "z_um", "y_um", "x_um", "effect_size",
                   "spot_mean_intensity", "detector"]


class TestSetGuardError(PermissionError):
    """A test gonad was asked for without ``allow_test`` (only ``bench_detector.py --final`` sets it)."""

    __test__ = False        # not a pytest test class, whatever its name says


# ----------------------------------------------------------------------------------- roles
def preregistered_role(tag: str) -> str | None:
    """The role the pre-registration gives ``tag`` (train, val or test), or None for a tag outside the
    N2 RAD-51 split (a synthetic or future gonad, whose role then comes from its manifest)."""
    from germquant.validate.benchmark_split import PREREGISTERED_ROLES

    return PREREGISTERED_ROLES.get(str(tag))


def gonad_role(tag: str, manifest_role: str | None = None) -> str:
    """The role of a gonad: the pre-registered one when the tag is in the split, else the manifest's.
    ValueError when both exist and disagree (a manifest that moves a gonad is not the registered split)."""
    pre = preregistered_role(tag)
    man = str(manifest_role or "").strip()
    if pre is not None and man and man != pre:
        raise ValueError(f"{tag}: the manifest says role {man!r} but the pre-registered split says {pre!r}")
    return pre or man


def entry_role(entry: Entry) -> str:
    """The role of a manifest entry (``role``, else ``split``), checked against the pre-registration."""
    return gonad_role(entry.tag, entry.extra.get("role") or entry.split)


def is_test_tag(tag: str, manifest_role: str | None = None) -> bool:
    """True when either the pre-registration or the manifest calls the gonad a test gonad."""
    return preregistered_role(tag) == "test" or str(manifest_role or "").strip() == "test"


def guard_test(tags, *, allow_test: bool = False, roles: dict | None = None) -> None:
    """Raise TestSetGuardError when any of ``tags`` is a test gonad and ``allow_test`` is false.
    ``roles`` maps a tag to its manifest role, for tags outside the pre-registered split."""
    roles = roles or {}
    test = [t for t in tags if is_test_tag(t, roles.get(t))]
    if test and not allow_test:
        raise TestSetGuardError(
            f"refused: {', '.join(test)} {'is a test gonad' if len(test) == 1 else 'are test gonads'} "
            "(docs/BENCHMARK_PROTOCOL.md). Test gonads are scored only by scripts/bench_detector.py --final "
            "with a frozen method and a --reason.")


# ------------------------------------------------------------------------------ small helpers
def canonical_json(obj) -> str:
    """Key-sorted compact JSON, the form whose sha256 identifies a parameter set."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str)


def params_sha256(params) -> str:
    """sha256 of the canonical JSON of a parameter mapping (``{}`` for none)."""
    return hashlib.sha256(canonical_json(params or {}).encode("utf-8")).hexdigest()


def file_sha256(path, chunk_bytes: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_bytes), b""):
            h.update(chunk)
    return h.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _jsonable(obj):
    from .report import jsonable

    return jsonable(obj)


def _write_json_atomic(path: Path, doc: dict) -> None:
    tmp = path.with_name(path.name + ".partial")
    tmp.write_text(json.dumps(_jsonable(doc), indent=1) + "\n", encoding="utf-8")
    os.replace(tmp, path)


# ------------------------------------------------------------------------------ cache writing
def write_npy_chunked(path, arr, z_chunk: int = DEFAULT_Z_CHUNK) -> str:
    """Write a (Z, ...) array to ``path`` as .npy through ``numpy.lib.format.open_memmap``, ``z_chunk``
    planes at a time, under a temporary name that is renamed into place only when complete. Returns the
    sha256 of the array data in C order (the bytes of the planes, not of the file with its header)."""
    path = Path(path)
    arr = np.asarray(arr)
    if arr.ndim < 1:
        raise ValueError("write_npy_chunked needs an array with at least one axis")
    z_chunk = max(1, int(z_chunk))
    tmp = path.with_name(path.stem + ".partial.npy")
    h = hashlib.sha256()
    mm = np.lib.format.open_memmap(tmp, mode="w+", dtype=arr.dtype, shape=arr.shape)
    try:
        for z0 in range(0, arr.shape[0], z_chunk):
            block = np.ascontiguousarray(arr[z0:z0 + z_chunk])
            mm[z0:z0 + z_chunk] = block
            h.update(block)
        mm.flush()
    finally:
        del mm
    os.replace(tmp, path)
    return h.hexdigest()


def array_sha256(arr, z_chunk: int = DEFAULT_Z_CHUNK) -> str:
    """sha256 of an array's data in C order, read ``z_chunk`` planes at a time (works on a memmap)."""
    h = hashlib.sha256()
    for z0 in range(0, arr.shape[0], max(1, int(z_chunk))):
        h.update(np.ascontiguousarray(arr[z0:z0 + z_chunk]))
    return h.hexdigest()


def label_voxel_counts(labels) -> np.ndarray:
    """Voxels per label value (index = label), counted plane by plane so a 2600 x 2600 x 71 int32
    image never needs a whole-volume int64 copy."""
    top = int(labels.max()) if labels.size else 0
    if labels.size and int(labels.min()) < 0:
        raise ValueError("label images must not hold negative labels")
    counts = np.zeros(top + 1, np.int64)
    for z in range(labels.shape[0]):
        counts += np.bincount(np.asarray(labels[z]).ravel(), minlength=top + 1)
    return counts


def golden_paths(results_dir, stem: str) -> dict[str, Path]:
    """The golden run's files under ``results_dir`` for an image id ``stem``."""
    rdir = Path(results_dir)
    return {k: rdir / pat.format(stem=stem) for k, pat in GOLDEN_FILES.items()}


def _entry_stem(entry: Entry) -> str:
    if entry.stem:
        return entry.stem
    if entry.nd2 is not None:
        return Path(entry.nd2).stem
    raise ValueError(f"{entry.tag}: the manifest entry names neither a stem nor an .nd2")


def role_indices_from_summary(summary: pd.DataFrame) -> dict:
    """{role: channel index or None} from the ``channel_map_json`` column of an image summary."""
    if summary.empty or "channel_map_json" not in summary.columns:
        raise ValueError("the image summary has no channel_map_json column")
    raw = json.loads(str(summary["channel_map_json"].iloc[0]))
    return {str(k): (None if v is None else int(v)) for k, v in raw.items()}


def _marker_from_tables(spots_path: Path) -> str | None:
    try:
        col = pd.read_csv(spots_path, usecols=["marker"])["marker"].dropna().astype(str).unique()
    except (ValueError, KeyError):
        return None
    return str(col[0]) if len(col) == 1 else None


def _marker_from_channel_map(run_manifest: dict) -> str | None:
    """The foci marker of the channel map the golden run's config names, resolved against this
    repository (the shipped channel maps are never edited, so the file is the one the run read)."""
    rel = ((run_manifest.get("config") or {}).get("io") or {}).get("channel_map")
    if not rel:
        return None
    repo = Path(__file__).resolve().parents[3]
    path = Path(rel) if Path(rel).is_absolute() else repo / rel
    if not path.is_file():
        return None
    from germquant.io.channel_map import ChannelMap

    return ChannelMap.from_yaml(path).marker("foci")


def _split_gonad(split_record: dict | None, tag: str) -> dict | None:
    if not split_record:
        return None
    for g in split_record.get("gonads", []):
        if g.get("tag") == tag:
            return g
    return None


def cache_gonad(entry: Entry, cache_root=DEFAULT_CACHE_ROOT, *, split_record: dict | None = None,
                reader=None, z_chunk: int = DEFAULT_Z_CHUNK, verify_sha256: bool = False,
                overwrite: bool = False) -> dict:
    """Write ``<cache_root>/<tag>/`` for one manifest entry and return a status record
    (``status`` is "written", "exists" or "skipped", with ``reason`` for a skip).

    ``reader`` reads the image (default ``germquant.io.read.read_image``, the reader the pipeline uses);
    ``split_record`` is the pre-registered split (tests/golden/rad51_split.json) that supplies the
    .nd2 path, size and sha256 for the tag. A results folder without the pipeline outputs is skipped,
    not fatal. ``verify_sha256`` hashes the .nd2 and refuses it when it is not the registered file.
    """
    tag = entry.tag
    out = Path(cache_root) / tag
    meta_path = out / CACHE_FILES["meta"]
    if meta_path.is_file() and not overwrite:
        return {"tag": tag, "status": "exists", "path": str(out)}
    role = gonad_role(tag, entry.extra.get("role") or entry.split)    # a role conflict fails before any write
    stem = _entry_stem(entry)
    paths = golden_paths(entry.results_dir, stem)
    missing = [k for k in REQUIRED_GOLDEN if not paths[k].is_file()]
    if missing:
        return {"tag": tag, "status": "skipped", "path": str(out),
                "reason": f"no pipeline outputs in {entry.results_dir} (missing {', '.join(missing)})"}

    run_manifest = json.loads(paths["run_manifest"].read_text(encoding="utf-8"))
    spots_cfg = copy.deepcopy((run_manifest.get("config") or {}).get("spots") or {})
    summary = pd.read_csv(paths["image_summary"])
    roles = role_indices_from_summary(summary)
    foci_idx, dna_idx = roles.get("foci"), roles.get("dna")
    if foci_idx is None:
        return {"tag": tag, "status": "skipped", "path": str(out),
                "reason": "the golden run resolved no foci channel (channel_map_json)"}
    nuclei_all = pd.read_csv(paths["nuclei"])
    marker = _marker_from_tables(paths["spots"])
    marker_source = "golden spots table"
    cm_marker = _marker_from_channel_map(run_manifest)
    if marker is None:
        marker, marker_source = cm_marker, "channel map of the golden config"
    elif cm_marker is not None and cm_marker != marker:
        raise ValueError(f"{tag}: the golden spots table says marker {marker!r} but the channel map says "
                         f"{cm_marker!r}")
    if marker is None:
        raise ValueError(f"{tag}: no marker name in the golden spots table or the channel map")

    # the image, from the registered .nd2
    g = _split_gonad(split_record, tag)
    nd2_rec = (g or {}).get("nd2") or {}
    nd2_path = Path(nd2_rec.get("path") or entry.nd2 or "")
    if not nd2_path.is_file():
        return {"tag": tag, "status": "skipped", "path": str(out), "reason": f"image not found: {nd2_path}"}
    size = nd2_path.stat().st_size
    if nd2_rec.get("size_bytes") is not None and size != int(nd2_rec["size_bytes"]):
        raise ValueError(f"{tag}: {nd2_path} is {size} bytes, the split registered {nd2_rec['size_bytes']}")
    sha_verified = False
    if verify_sha256 and nd2_rec.get("sha256"):
        got = file_sha256(nd2_path)
        if got != nd2_rec["sha256"]:
            raise ValueError(f"{tag}: {nd2_path} sha256 {got} is not the registered {nd2_rec['sha256']}")
        sha_verified = True

    import tifffile

    labels = np.asarray(tifffile.imread(str(paths["nuclei_labels"])))
    if labels.ndim == 2:
        labels = labels[None]
    if labels.dtype.kind not in "iu":
        raise ValueError(f"{tag}: the golden label image is {labels.dtype}, not an integer image")
    if labels.dtype != np.int32:
        if labels.size and (int(labels.min()) < 0 or int(labels.max()) > np.iinfo(np.int32).max):
            raise ValueError(f"{tag}: label values do not fit int32")
        labels = labels.astype(np.int32)
    voxels_per_label = label_voxel_counts(labels)
    label_ids = np.flatnonzero(voxels_per_label[1:]) + 1
    table_ids = set(nuclei_all["nucleus_id"].astype(int).tolist())
    if set(label_ids.tolist()) != table_ids:
        raise ValueError(f"{tag}: the label image holds {len(label_ids)} labels, the nuclei table "
                         f"{len(table_ids)} ids, and the two sets differ")

    if reader is None:
        from germquant.io.read import read_image as reader
    t_read = time.time()
    stack = reader(nd2_path)
    read_s = time.time() - t_read
    data = stack.data
    if data.ndim != 4:
        raise ValueError(f"{tag}: the reader returned {data.shape}, expected (C, Z, Y, X)")
    if tuple(data.shape[1:]) != tuple(labels.shape):
        raise ValueError(f"{tag}: image {tuple(data.shape[1:])} and label image {labels.shape} differ in shape")
    spacing = tuple(float(s) for s in stack.spacing)
    tvox = [c for c in ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um") if c in nuclei_all.columns]
    if len(tvox) == 3:
        tv = tuple(float(nuclei_all[c].iloc[0]) for c in tvox)
        if not np.allclose(tv, spacing, rtol=1e-9, atol=0):
            raise ValueError(f"{tag}: reader voxel size {spacing} differs from the golden tables' {tv}")
    if str(data.dtype) != "uint16":
        log.warning("%s: the reader returned %s, not uint16; the cache keeps the reader's dtype so the "
                    "detector sees the same pixels", tag, data.dtype)

    out.mkdir(parents=True, exist_ok=True)
    if meta_path.exists():
        meta_path.unlink()          # an overwritten cache must never look complete half way through
    t_write = time.time()
    arrays = {}
    for name, idx in (("foci", foci_idx), ("dna", dna_idx)):
        if idx is None:
            continue
        sha = write_npy_chunked(out / CACHE_FILES[name], data[idx], z_chunk)
        arrays[name] = {"file": CACHE_FILES[name], "dtype": str(data.dtype), "shape": list(data.shape[1:]),
                        "channel_index": int(idx), "data_sha256": sha}
    channel_names = [str(c) for c in (getattr(stack, "channel_names", None) or [])]
    del stack, data
    sha = write_npy_chunked(out / CACHE_FILES["labels"], labels, z_chunk)
    arrays["labels"] = {"file": CACHE_FILES["labels"], "dtype": "int32", "shape": list(labels.shape),
                        "data_sha256": sha}
    write_s = time.time() - t_write

    keep = [c for c in NUCLEI_COLUMNS if c in nuclei_all.columns]
    nuc = nuclei_all[keep].copy().sort_values("nucleus_id", kind="stable").reset_index(drop=True)
    nuc["in_germline"] = bs.germline_mask(nuc) if "in_germline" in nuc.columns else True
    nuc.to_csv(out / CACHE_FILES["nuclei"], index=False)
    germ_ids = sorted(int(v) for v in nuc.loc[nuc["in_germline"], "nucleus_id"])

    from germquant import provenance

    golden_files = {}
    for k, p in paths.items():
        if p.is_file():
            golden_files[k] = {"path": str(p).replace("\\", "/"), "size_bytes": p.stat().st_size}
    meta = {
        "format_version": CACHE_FORMAT,
        "tag": tag,
        "role": role,
        "condition": entry.extra.get("condition"),
        "sex": entry.extra.get("sex"),
        "image_id": str(summary["image_id"].iloc[0]) if "image_id" in summary.columns else stem,
        "stem": stem,
        "spacing_zyx_um": list(spacing),
        "shape_zyx": list(labels.shape),
        "role_indices": roles,
        "channel_names": channel_names,
        "marker": marker,
        "marker_source": marker_source,
        "spots_config": spots_cfg,
        "nd2": {"path": str(nd2_path).replace("\\", "/"), "size_bytes": size,
                "sha256": nd2_rec.get("sha256"), "sha256_source": "split record" if nd2_rec.get("sha256") else None,
                "sha256_verified": sha_verified},
        "labels": {
            "source": str(paths["nuclei_labels"]).replace("\\", "/"),
            "content": "every segmented nucleus (the label image the spots stage received; germline "
                       "and off-gonad alike)",
            "n_labels": int(len(label_ids)), "max_label": int(label_ids.max()) if len(label_ids) else 0,
        },
        "n_nuclei": int(len(nuc)),
        "n_germline": len(germ_ids),
        "germline_label_ids": germ_ids,
        "arrays": arrays,
        "golden": {
            "results_dir": str(entry.results_dir).replace("\\", "/"),
            "git_sha": run_manifest.get("git_sha"),
            "config_hash": run_manifest.get("config_hash"),
            "run_timestamp": run_manifest.get("run_timestamp"),
            "tool_versions": run_manifest.get("tool_versions"),
            "n_spots": int(len(pd.read_csv(paths["spots"], usecols=["spot_id"]))),
            "files": golden_files,
        },
        "created_utc": _utc_now(),
        "germquant_git_sha": provenance.git_sha(),
        "z_chunk": int(z_chunk),
        "timing_s": {"read": round(read_s, 1), "write": round(write_s, 1)},
    }
    _write_json_atomic(meta_path, meta)
    nbytes = sum((out / f).stat().st_size for f in CACHE_FILES.values() if (out / f).is_file())
    return {"tag": tag, "status": "written", "path": str(out), "bytes": int(nbytes)}


def cache_size_bytes(cache_root=DEFAULT_CACHE_ROOT) -> dict:
    """Bytes on disk per cached gonad (complete caches only) and in total."""
    root = Path(cache_root)
    per = {}
    if root.is_dir():
        for d in sorted(p for p in root.iterdir() if p.is_dir()):
            if (d / CACHE_FILES["meta"]).is_file():
                per[d.name] = int(sum(f.stat().st_size for f in d.iterdir() if f.is_file()))
    return {"per_tag": per, "total": int(sum(per.values()))}


# ------------------------------------------------------------------------------ cache reading
@dataclass
class GonadCache:
    tag: str
    path: Path
    meta: dict
    nuclei: pd.DataFrame
    foci: np.ndarray | None = None
    dna: np.ndarray | None = None
    labels: np.ndarray | None = None

    @property
    def spacing(self) -> tuple[float, float, float]:
        s = self.meta["spacing_zyx_um"]
        return (float(s[0]), float(s[1]), float(s[2]))

    @property
    def spots_config(self) -> dict:
        return copy.deepcopy(self.meta.get("spots_config") or {})

    @property
    def marker(self) -> str:
        return str(self.meta["marker"])

    @property
    def germline_ids(self) -> list[int]:
        return [int(v) for v in self.meta.get("germline_label_ids", [])]

    @property
    def role(self) -> str:
        return str(self.meta.get("role") or "")

    def array(self, name: str, *, mmap: bool = True) -> np.ndarray:
        """One cached array, loaded now (``mmap`` maps the file instead of reading it)."""
        return _load_array(self.path, self.meta, name, mmap=mmap)


def _load_array(path: Path, meta: dict, name: str, *, mmap: bool) -> np.ndarray:
    rec = (meta.get("arrays") or {}).get(name)
    if rec is None:
        raise FileNotFoundError(f"{meta.get('tag')}: the cache holds no {name} array")
    arr = np.load(path / rec["file"], mmap_mode="r" if mmap else None, allow_pickle=False)
    if list(arr.shape) != list(rec["shape"]) or str(arr.dtype) != rec["dtype"]:
        raise ValueError(f"{meta.get('tag')}: {rec['file']} is {arr.dtype} {arr.shape}, meta.json says "
                         f"{rec['dtype']} {rec['shape']}")
    return arr


def load_cache(tag: str, cache_root=None, *, arrays=ARRAYS, mmap: bool = True,
               verify: bool = False) -> GonadCache:
    """Read ``<cache_root>/<tag>/``. ``arrays`` names the arrays to open (any of foci, dna, labels);
    ``mmap`` maps them read-only instead of reading them into memory; ``verify`` recomputes the data
    sha256 of every opened array and compares it with meta.json. FileNotFoundError when the cache is
    missing or unfinished (no meta.json)."""
    path = Path(cache_root or DEFAULT_CACHE_ROOT) / tag
    meta_path = path / CACHE_FILES["meta"]
    if not meta_path.is_file():
        raise FileNotFoundError(f"no complete cache for {tag} under {path.parent} (scripts/bench_cache.py "
                                "writes it; meta.json is written last)")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if meta.get("format_version") != CACHE_FORMAT:
        raise ValueError(f"{tag}: cache format {meta.get('format_version')!r}, this code reads {CACHE_FORMAT}")
    nuclei = pd.read_csv(path / CACHE_FILES["nuclei"])
    if "in_germline" in nuclei.columns:
        nuclei["in_germline"] = bs.germline_mask(nuclei)
    out = GonadCache(tag=tag, path=path, meta=meta, nuclei=nuclei)
    for name in arrays:
        if name not in ARRAYS:
            raise ValueError(f"unknown cache array {name!r}; the cache holds {ARRAYS}")
        if name not in (meta.get("arrays") or {}):
            continue
        arr = _load_array(path, meta, name, mmap=mmap)
        if verify and array_sha256(arr, meta.get("z_chunk", DEFAULT_Z_CHUNK)) != meta["arrays"][name]["data_sha256"]:
            raise ValueError(f"{tag}: {name} data sha256 differs from meta.json")
        setattr(out, name, arr)
    return out


# --------------------------------------------------------------------------------- running
def effective_section(spots_config: dict | None, params: dict | None) -> dict:
    """The ``spots:`` section a run uses: the cached golden section with ``params`` laid over it. A
    nested mapping (``log:``) merges one level deep over the section's own, the rule
    ``pipeline.spot_backend_section`` applies to a spot instance; every other key replaces."""
    section = copy.deepcopy(dict(spots_config or {}))
    params = copy.deepcopy(dict(params or {}))
    if "backend" in params:
        raise ValueError("pass the backend as its own argument (--backend), not inside params")
    for k, v in params.items():
        base = section.get(k)
        section[k] = {**base, **v} if isinstance(v, dict) and isinstance(base, dict) else v
    return section


def pipeline_spot_kwargs(section: dict) -> dict:
    """The SpotMAX keyword arguments the pipeline's spots stage builds from its ``spots:`` section, with
    the same expressions and code defaults (src/germquant/pipeline.py, ``_spots``)."""
    s = section
    return dict(
        spot_radius_um=float(s.get("spot_radius_um", 0.3)),
        gauss_sigma_um=float(s.get("gauss_sigma_um", 0.08)),
        thresholding_method=s.get("thresholding_method", "threshold_triangle"),
        effect_size_metric=s.get("effect_size_metric", "spot_vs_backgr_effect_size_glass"),
        effect_size_min=float(s.get("effect_size_min", 3.0)),
        merge_z_columns=bool(s.get("merge_z_columns", True)),
        z_merge_gap_um=float(s.get("z_merge_gap_um", 0.8)),
        z_merge_valley_frac=float(s.get("z_merge_valley_frac", 0.8)),
        max_spot_candidates=int(s.get("max_spot_candidates", 30000)),
    )


def _in_memory(arr):
    """A plain in-memory ndarray (the pipeline hands the detector arrays, not file maps)."""
    if arr is None:
        return None
    if isinstance(arr, np.memmap):
        return np.array(arr)
    return np.asarray(arr)


# channel roles a backend may ask for (its input_roles) that the cache holds, and the cached array of each
CACHED_ROLES = {"dna": "dna"}


def backend_cache_arrays(backend: str) -> tuple[str, ...]:
    """The cache arrays a run of ``backend`` reads: foci and labels, plus the cached channel of every
    extra role it asks for (net3d: dna). ValueError for an unknown backend, as detect_spots."""
    from germquant.spots.backends import backend_input_roles

    extra = tuple(CACHED_ROLES[r] for r in backend_input_roles(backend) if r in CACHED_ROLES)
    return ("foci", "labels") + extra


def run_backend(tag: str, backend: str = "spotmax", params: dict | None = None, *, cache_root=None,
                cache: GonadCache | None = None, allow_test: bool = False):
    """Run spot detection backend ``backend`` on the cached gonad ``tag`` exactly as the pipeline's
    spots stage calls ``detect_spots`` (same image, labels, spacing, marker, keyword parameters from the
    cached spots section with ``params`` overriding, ``backend`` and ``cfg_section``; for a backend that
    asks for further channel roles, ``channels`` with those the cache holds, as the stage passes those
    the image has). Returns (per_spot, per_nucleus), the two tables of ``detect_spots``. A test gonad
    is refused unless ``allow_test``."""
    from germquant.spots import detect_spots
    from germquant.spots.backends import backend_input_roles

    guard_test([tag], allow_test=allow_test, roles={tag: cache.role if cache is not None else None})
    roles = backend_input_roles(backend)
    if cache is None:
        cache = load_cache(tag, cache_root, arrays=backend_cache_arrays(backend), mmap=False)
    guard_test([tag], allow_test=allow_test, roles={tag: cache.role})
    foci = cache.foci if cache.foci is not None else cache.array("foci", mmap=False)
    labels = cache.labels if cache.labels is not None else cache.array("labels", mmap=False)
    section = effective_section(cache.spots_config, params)
    kw = pipeline_spot_kwargs(section)
    extra = {}
    if roles:
        chans = {}
        for r in roles:
            name = CACHED_ROLES.get(r)
            if name is None or name not in (cache.meta.get("arrays") or {}):
                continue                      # the role is not in the cache: like an image without it
            arr = getattr(cache, name, None)
            chans[r] = _in_memory(arr if arr is not None else cache.array(name, mmap=False))
        extra = {"channels": chans}
    return detect_spots(_in_memory(foci), _in_memory(labels), cache.spacing, marker=cache.marker, **kw,
                        backend=backend, cfg_section=section, **extra)


# --------------------------------------------------------------------------------- scoring
def _entries_by_tag(manifest) -> dict[str, Entry]:
    if manifest is None:
        manifest = DEFAULT_MANIFEST
    if isinstance(manifest, (str, Path)):
        entries = load_manifest(manifest)
    else:
        entries = list(manifest)
    return {e.tag: e for e in entries}


def spot_truth_usable(entry: Entry) -> bool:
    """False for a gonad whose spots were exported with Imaris region growing (manifest flag)."""
    return entry.extra.get("spot_truth_usable") is not False


def score(tag: str, per_spot: pd.DataFrame, per_nucleus: pd.DataFrame, *, manifest=None, cache_root=None,
          cache: GonadCache | None = None, truth=None, truth_labels=None, allow_test: bool = False,
          include_unusable: bool = False, tol_xy_um: float = bs.DEFAULT_TOL_XY_UM,
          tol_z_um: float = bs.DEFAULT_TOL_Z_UM, surface_tol_um: float = bs.DEFAULT_SURFACE_TOL_UM,
          overlap_frac: float = bs.DEFAULT_OVERLAP_FRAC, traced_rule: str = "auto"):
    """Score one run of a detector on a cached gonad against the truth files the manifest lists.

    Reuses ``germquant.benchmark.spots.benchmark_gonad`` with the detector's spots and per-nucleus
    counts in place of the golden ones and the cached label image as our nuclei: Hungarian matching in
    the ``tol_xy_um`` by ``tol_z_um`` cylinder (all expert spots, without Class B, and detector-only in
    the traced nuclei), per-nucleus agreement of the expert spots mapped into our labels with our
    counts over the germline nuclei and the traced subset. On top of that the per-nucleus agreement
    with the truth counted from Class A spots only (``class_a_*`` columns; the expert count changes,
    ours does not, as section 5.1 of the protocol states it). Returns (per-gonad row, per-nucleus
    table, matched pairs). Refuses a test gonad unless ``allow_test``, and a gonad whose spot truth
    is unusable unless ``include_unusable``."""
    entries = _entries_by_tag(manifest)
    if tag not in entries:
        raise KeyError(f"{tag} is not in the manifest")
    entry = entries[tag]
    role = entry_role(entry)
    guard_test([tag], allow_test=allow_test, roles={tag: role})
    if not spot_truth_usable(entry) and not include_unusable:
        raise ValueError(f"{tag}: the spot truth is unusable (exported with Imaris region growing)")
    if cache is None:
        cache = load_cache(tag, cache_root, arrays=("labels",), mmap=True)
    labels = cache.labels if cache.labels is not None else cache.array("labels", mmap=True)
    if truth is None:
        truth = load_truth(entry.truth_dir, entry.truth_stem)
        if truth is None:
            raise FileNotFoundError(f"{tag}: no truth files for stem {entry.truth_stem!r} in {entry.truth_dir}")
    if truth_labels is None and truth.has_labels:
        truth_labels = load_labels(truth.labels_path)

    nuclei = cache.nuclei.copy()
    counts = (per_nucleus.set_index("nucleus_id")["n_spots"] if len(per_nucleus)
              else pd.Series(dtype=float))
    nuclei["n_spots"] = nuclei["nucleus_id"].map(counts).fillna(0).astype(int)
    unknown = sorted(set(per_nucleus["nucleus_id"].astype(int)) - set(nuclei["nucleus_id"].astype(int))) \
        if len(per_nucleus) else []
    spacing = cache.spacing
    outputs = PipelineOutputs(cache.meta.get("image_id", tag), nuclei, per_spot, None, spacing, "bench cache")
    condition = str(entry.extra.get("condition") or cache.meta.get("condition") or "")
    row, per_nuc, matches = bs.benchmark_gonad(
        outputs, truth, labels, spacing, tag=tag, split=condition, truth_labels=truth_labels,
        tol_xy_um=tol_xy_um, tol_z_um=tol_z_um, surface_tol_um=surface_tol_um, overlap_frac=overlap_frac,
        traced_rule=traced_rule)
    row["voxel_source"] = "bench cache"
    row["role"] = role
    row["n_unknown_nucleus_ids"] = len(unknown)

    # Class A only: the expert count from Class A spots, ours unchanged (protocol section 5.1)
    nan = float("nan")
    if truth.has_spot_class and len(truth.spots):
        a = ~bs.is_class_b(truth.spots["spot_class"])
        tz, ty, tx = bs._spot_voxels(truth.spots, spacing)
        lab_a = bs.label_at_points(labels, tz[a], ty[a], tx[a])
        counts_a = pd.Series(lab_a[lab_a > 0]).value_counts()
        per_nuc["expert_n_class_a"] = per_nuc["nucleus_id"].map(counts_a).fillna(0).astype(int)
        row["n_expert_spots_class_a"] = int(a.sum())
        m_all = agreement_metrics(per_nuc["expert_n_class_a"], per_nuc["ours_n"])
        tr = per_nuc[per_nuc["traced"]]
        m_tr = agreement_metrics(tr["expert_n_class_a"], tr["ours_n"])
        for k in ("mean_expert", "pearson_r", "ccc", "bias"):
            row[f"class_a_{k}_all"] = m_all[k]
            row[f"class_a_{k}_traced"] = m_tr[k]
    else:
        per_nuc["expert_n_class_a"] = np.nan
        row["n_expert_spots_class_a"] = 0
        for k in ("mean_expert", "pearson_r", "ccc", "bias"):
            row[f"class_a_{k}_all"] = nan
            row[f"class_a_{k}_traced"] = nan
    per_nuc.insert(2, "role", role)
    return row, per_nuc, matches


def pool_scores(rows: pd.DataFrame, per_nuc: pd.DataFrame, *, n_boot: int = N_BOOT,
                seed: int = BOOT_SEED) -> dict:
    """Pooled scores over the scored gonads with cluster-bootstrap intervals by gonad: the blocks of
    ``germquant.benchmark.spots.pool`` (by condition in ``by_split`` when both are present) plus the
    per-nucleus agreement with the Class A truth counts (``per_nucleus_all_class_a``,
    ``per_nucleus_traced_class_a``)."""
    out = bs.pool(rows, per_nuc, n_boot=n_boot, seed=seed)

    def class_a_blocks(r: pd.DataFrame, pn: pd.DataFrame) -> dict:
        if "expert_n_class_a" not in pn.columns or not pn["expert_n_class_a"].notna().any():
            return {}
        pa = pn[pn["expert_n_class_a"].notna()].assign(expert_n=lambda d: d["expert_n_class_a"])
        tags = [t for t in r["tag"] if t in set(pa["tag"])]
        return {"per_nucleus_all_class_a": bs._agreement_block(pa, tags, n_boot=n_boot, seed=seed),
                "per_nucleus_traced_class_a": bs._agreement_block(pa[pa["traced"]], tags, n_boot=n_boot,
                                                                  seed=seed)}

    out.update(class_a_blocks(rows, per_nuc))
    for s, block in (out.get("by_split") or {}).items():
        block.update(class_a_blocks(rows[rows["split"].astype(str) == s],
                                    per_nuc[per_nuc["split"].astype(str) == s]))
    return out


# ----------------------------------------------------------------------- golden comparison
def read_golden_spots(results_dir, stem: str) -> pd.DataFrame:
    """The golden per-spot table, from the parquet file when it exists (exact binary floats), else from
    the CSV parsed with round-trip float precision (pandas wrote the shortest exact representation)."""
    p = golden_paths(results_dir, stem)
    if p["spots_parquet"].is_file():
        try:
            return pd.read_parquet(p["spots_parquet"])
        except Exception as e:  # noqa: BLE001 - pyarrow missing: the CSV is exact too
            log.warning("could not read %s (%s); using the CSV", p["spots_parquet"], e)
    return pd.read_csv(p["spots"], float_precision="round_trip")


def compare_with_golden(per_spot: pd.DataFrame, per_nucleus: pd.DataFrame, golden_spots: pd.DataFrame,
                        golden_nuclei: pd.DataFrame | None = None) -> dict:
    """Column by column comparison of a detector run with the golden spots table (row order included)
    and of the per-nucleus counts with the golden ``n_spots``. ``identical`` is true when every compared
    column is exactly equal (floats bit for bit, NaN equal to NaN) and the counts agree."""
    res: dict = {"n_ours": int(len(per_spot)), "n_golden": int(len(golden_spots)), "columns": {}}
    same_len = len(per_spot) == len(golden_spots)
    ok = same_len
    for c in COMPARE_COLUMNS:
        if c not in per_spot.columns or c not in golden_spots.columns:
            res["columns"][c] = {"compared": False, "reason": "column missing"}
            ok = ok and c not in golden_spots.columns
            continue
        if not same_len:
            res["columns"][c] = {"compared": False, "reason": "row counts differ"}
            continue
        a = per_spot[c].to_numpy()
        b = golden_spots[c].to_numpy()
        if a.dtype.kind in "fc" or b.dtype.kind in "fc":
            af, bf = a.astype(float), b.astype(float)
            eq = (af == bf) | (np.isnan(af) & np.isnan(bf))
            with np.errstate(invalid="ignore"):          # inf - inf where both sides are infinite
                diff = np.abs(af - bf)
            diff = diff[np.isfinite(diff)]
            rec = {"compared": True, "n_different": int((~eq).sum()),
                   "max_abs_diff": float(diff.max()) if diff.size else 0.0}
        else:
            eq = np.asarray([str(x) == str(y) for x, y in zip(a, b)]) if a.dtype.kind == "O" or b.dtype.kind == "O" \
                else (a == b)
            rec = {"compared": True, "n_different": int((~np.asarray(eq)).sum())}
        res["columns"][c] = rec
        ok = ok and rec["n_different"] == 0
    if golden_nuclei is not None and "n_spots" in golden_nuclei.columns:
        mine = per_nucleus.set_index("nucleus_id")["n_spots"] if len(per_nucleus) else pd.Series(dtype=int)
        gold = golden_nuclei.set_index("nucleus_id")["n_spots"].fillna(0).astype(int)
        ids = sorted(set(gold.index.astype(int)) | set(mine.index.astype(int)))
        m = np.array([int(mine.get(i, 0)) for i in ids])
        g = np.array([int(gold.get(i, 0)) for i in ids])
        res["nuclei"] = {"n_nuclei": len(ids), "n_count_different": int((m != g).sum()),
                         "sum_ours": int(m.sum()), "sum_golden": int(g.sum())}
        ok = ok and res["nuclei"]["n_count_different"] == 0
    res["identical"] = bool(ok)
    return res


# ------------------------------------------------------------------------------- trial log
POOLED_KEYS = ("per_nucleus_all", "per_nucleus_traced", "per_nucleus_all_class_a", "per_nucleus_traced_class_a",
               "spots_all", "spots_without_class_b", "spots_detector_only")
PER_GONAD_KEYS = ("n_our_spots", "n_expert_spots", "ccc_all", "bias_all", "ccc_traced", "bias_traced",
                  "class_a_ccc_all", "class_a_ccc_traced", "precision", "recall", "f1", "no_b_f1", "det_f1")


def trial_record(*, backend: str, params: dict, split: str, tags: list[str], pooled: dict,
                 per_gonad: pd.DataFrame, skipped: list, out_dir, manifest, cache_root, n_boot: int,
                 seed: int, elapsed_s: float, final: bool = False, reason: str | None = None) -> dict:
    """One line of trials.jsonl (or final_runs.jsonl): when, which code, which detector and parameters,
    which gonads, and the pooled metrics with their intervals."""
    from germquant import provenance

    versions = provenance.tool_versions()
    rec = {
        "utc": _utc_now(),
        "git_sha": provenance.git_sha(),
        "backend": backend,
        "params": params or {},
        "params_sha256": params_sha256(params),
        "split": split,
        "tags": list(tags),
        "skipped": skipped,
        "final": bool(final),
        "pooled": {k: pooled[k] for k in POOLED_KEYS if k in pooled},
        "pooled_by_condition": {s: {k: b[k] for k in POOLED_KEYS if k in b}
                                for s, b in (pooled.get("by_split") or {}).items()},
        "per_gonad": {str(r["tag"]): {k: r[k] for k in PER_GONAD_KEYS if k in r}
                      for _, r in per_gonad.iterrows()},
        "bootstrap": {"n_boot": int(n_boot), "seed": int(seed)},
        "out_dir": str(out_dir),
        "manifest": str(manifest),
        "cache_root": str(cache_root),
        "versions": {k: versions.get(k) for k in ("germquant", "spotmax", "numpy", "scipy", "scikit-image",
                                                  "pandas")},
        "elapsed_s": round(float(elapsed_s), 1),
    }
    if reason is not None:
        rec["reason"] = reason
    return rec


def append_jsonl(path, record: dict) -> None:
    """Append one JSON line (NaN written as null) in a single write."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(_jsonable(record), sort_keys=False) + "\n"
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(line)


def read_jsonl(path) -> list[dict]:
    path = Path(path)
    if not path.is_file():
        return []
    return [json.loads(s) for s in path.read_text(encoding="utf-8").splitlines() if s.strip()]


def n_final_runs(path, backend: str) -> int:
    """How many test-set scorings final_runs.jsonl already records for a method family (the backend)."""
    return sum(1 for r in read_jsonl(path) if r.get("backend") == backend)


# ---------------------------------------------------------------------------- one gonad
def bench_one(tag: str, backend: str, params: dict | None, *, cache_root=None, manifest=None,
              allow_test: bool = False, score_kwargs: dict | None = None) -> dict:
    """Run and score one gonad (the unit of work of scripts/bench_detector.py; importable, so a process
    pool can run it). Returns the per-gonad row, the per-nucleus table, the matched pairs and the
    detector's own per-spot table."""
    t0 = time.time()
    guard_test([tag], allow_test=allow_test)
    cache = load_cache(tag, cache_root, arrays=backend_cache_arrays(backend), mmap=False)
    per_spot, per_nuc = run_backend(tag, backend, params, cache=cache, allow_test=allow_test)
    t_detect = time.time() - t0
    cache.foci = None
    cache.dna = None
    t1 = time.time()
    row, pn, matches = score(tag, per_spot, per_nuc, manifest=manifest, cache=cache, allow_test=allow_test,
                             **(score_kwargs or {}))
    row["backend"] = backend
    row["detect_s"] = round(t_detect, 1)
    row["score_s"] = round(time.time() - t1, 1)
    row["n_spots_reconcile"] = bool(int(per_nuc["n_spots"].sum()) == len(per_spot)) if len(per_nuc) else len(per_spot) == 0
    return {"tag": tag, "row": row, "per_nucleus": pn, "matches": matches, "per_spot": per_spot}


def finite_or_none(x):
    """A float for JSON and printing, None for NaN."""
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None
