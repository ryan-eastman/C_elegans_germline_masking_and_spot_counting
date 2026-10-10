"""The pre-registered RAD-51 benchmark split and the manifests built from it.

docs/BENCHMARK_PROTOCOL.md is the pre-registration; this module is its machine-readable side.

* ``PREREGISTERED_ROLES`` fixes which of the 12 N2 RAD-51 gonads of 2025-11-05 (the ``n2rad51_*``
  tags of tests/golden/real_manifest.json) are train, val and test. ``tests/golden/rad51_split.json``
  is the record built from it by scripts/make_rad51_split.py: for every gonad the .nd2 stack and the
  Imaris .ims project with size and sha256, condition, sex, the Imaris object counts and whether the
  spots were exported with Imaris region growing on (read from the creation parameters stored in the
  .ims). ``load_split`` reads and validates that file; ``check_split_files`` confirms the files are
  still the ones that were registered.
* ``build_manifest`` joins the split with a pipeline results root (one folder per tag, as
  scripts/make_goldens.py writes them) and a truth root (the ``<stem>__imaris_*`` files of the Imaris
  truth export) and records, per gonad, every path and whether it exists. scripts/benchmark_manifest.py
  writes it to benchmarks/<name>/manifest.json; the scoring step reads that manifest.

Nothing here touches the pipeline: the module reads files and writes JSON. The Imaris reading is kept
separate from validate/imaris_ims.py on purpose, so the truth-export work can change that module
without touching the split record.
"""
from __future__ import annotations

import hashlib
import json
import re
from datetime import date, datetime, timezone
from pathlib import Path

import numpy as np

ROLES = ("train", "val", "test")

# The date the split was registered (docs/BENCHMARK_PROTOCOL.md, changelog entry 1).
REGISTERED_ON = "2026-10-09"

# The split. It never changes; see the protocol for why each gonad sits where it does.
PREREGISTERED_ROLES = {
    "n2rad51_noHS_herm_000": "train",
    "n2rad51_noHS_herm_001": "val",
    "n2rad51_noHS_herm_002": "test",
    "n2rad51_noHS_male_000": "train",
    "n2rad51_noHS_male_001": "train",
    "n2rad51_noHS_male_002": "test",
    "n2rad51_HS_herm_000": "train",
    "n2rad51_HS_herm_001": "val",
    "n2rad51_HS_herm_002": "test",
    "n2rad51_HS_male_000": "train",
    "n2rad51_HS_male_001": "test",
    "n2rad51_HS_male_002": "train",
}

TAG_RE = re.compile(r"^n2rad51_(?P<condition>noHS|HS)_(?P<sex>herm|male)_(?P<replicate>\d{3})$")

# Where the lab's Imaris projects sit under the ims root, per condition.
IMS_SUBDIR = {"noHS": "noHS_ims2", "HS": "HS_ims2"}

# Pipeline output files the benchmark reads, under <results root>/<tag>/, keyed by a short name.
# {stem} is the image id (the .nd2 file name without extension, spaces included).
OUTPUT_FILES = {
    "nuclei": "{stem}__nuclei.csv",
    "nuclei_labels": "{stem}__nuclei_labels.tif",
    "spots": "{stem}__spots.csv",
    "image_summary": "{stem}__image_summary.csv",
    "run_manifest": "run_manifest.json",
}
REQUIRED_OUTPUTS = ("nuclei", "nuclei_labels", "spots")

# Truth files of the Imaris export, under <truth root>, keyed the same way. The agreed contract is the
# name shape <stem>__imaris_<kind>.<ext>; the kinds and extensions below are the expected ones and can
# be overridden (build_manifest(truth_files=...), or --truth-file on the command line). Whatever else
# matches <stem>__imaris_* is listed under "found" so nothing the export wrote goes unnoticed.
TRUTH_FILES = {
    "spots": "{stem}__imaris_spots.csv",
    "surfaces": "{stem}__imaris_surfaces.csv",
    "masks": "{stem}__imaris_masks.tif",
    "meta": "{stem}__imaris_meta.json",
}
REQUIRED_TRUTH = ("spots", "surfaces", "masks")

NOTES = {
    "n2rad51_noHS_male_002": [
        "Floods the spot detector: SpotMAX with the deployed config returns far more candidates than the "
        "pipeline can process, so pipeline outputs are unavailable until a performance fix. The gonad is a "
        "legitimate test gonad on the truth side: it stays in the test set and is reported as missing, "
        "never dropped.",
    ],
}
REGION_GROWING_NOTE = (
    "Spots exported with Imaris region growing on (mEnableRegionGrowing is true in the Spots creation "
    "parameters). The spot truth of this gonad is excluded from spot scoring until it is re-exported "
    "without region growing; the surfaces are unaffected and remain nucleus truth.")
REGION_GROWING_UNKNOWN_NOTE = (
    "region_growing is null: the Spots creation parameters could not be read from the .ims, so whether "
    "region growing was on is unknown.")

_POINTS_CP = "Scene8/Content/Points0/CreationParameters"
_SPOTS = "Scene8/Content/Points0/Spot"
_SURF = "Scene8/Content/MegaSurfaces0/SurfaceModelInfo"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


# ----------------------------------------------------------------------------- files and hashes


def file_sha256(path, chunk_bytes: int = 1 << 26) -> str:
    """Streaming sha256 hex digest of a file (64 MiB chunks by default, so a 3 GB stack is hashed
    without being loaded into memory)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_bytes), b""):
            h.update(chunk)
    return h.hexdigest()


def posix(path) -> str:
    """Forward-slash form of a path, the spelling used throughout the manifests."""
    return Path(path).as_posix()


def file_record(path, sha256: str | None = None) -> dict:
    """{path, size_bytes, sha256} for a file; the digest is computed unless given."""
    p = Path(path)
    return {"path": posix(p), "size_bytes": p.stat().st_size, "sha256": sha256 or file_sha256(p)}


def load_hash_cache(path) -> dict:
    """Digests from an earlier run, keyed by posix path: either a plain {path: {size_bytes, sha256}}
    mapping or an existing split file (its nd2 and ims records are used)."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(data, dict) and "gonads" in data:
        cache = {}
        for g in data["gonads"]:
            for key in ("nd2", "ims"):
                rec = g.get(key) or {}
                if "path" in rec and "sha256" in rec:
                    cache[posix(rec["path"])] = {"size_bytes": rec.get("size_bytes"), "sha256": rec["sha256"]}
        return cache
    return {posix(k): v for k, v in data.items()}


def _record_with_cache(path, cache: dict | None) -> dict:
    """file_record, reusing a cached digest when the cached size still matches the file."""
    p = Path(path)
    hit = (cache or {}).get(posix(p))
    if hit and hit.get("size_bytes") == p.stat().st_size and hit.get("sha256"):
        return file_record(p, sha256=hit["sha256"])
    return file_record(p)


# ----------------------------------------------------------------------------- the Imaris project


def tag_parts(tag: str) -> dict:
    """condition (noHS|HS), sex (herm|male) and replicate (nnn) of an n2rad51 tag."""
    m = TAG_RE.match(tag)
    if not m:
        raise ValueError(f"{tag!r} is not an n2rad51 tag (n2rad51_<noHS|HS>_<herm|male>_<nnn>)")
    return m.groupdict()


def _attr_text(group, key: str) -> str | None:
    """Imaris stores attributes as arrays of single bytes; return them as text."""
    v = group.attrs.get(key)
    if v is None:
        return None
    if isinstance(v, np.ndarray) and v.dtype.kind == "S":
        return v.tobytes().decode("latin1").strip("\x00")
    if isinstance(v, bytes):
        return v.decode("latin1").strip("\x00")
    return str(v)


def _text_dataset(f, key: str) -> str | None:
    if key not in f:
        return None
    raw = f[key][()]
    if isinstance(raw, np.ndarray):
        raw = raw.ravel()[0] if raw.size else b""
    if isinstance(raw, (bytes, np.bytes_)):
        return bytes(raw).decode("latin1").strip("\x00")
    return str(raw)


def _xml_flag(text: str | None, name: str) -> bool | None:
    m = re.search(rf'\b{name}="(true|false)"', text or "")
    return None if m is None else m.group(1) == "true"


def read_ims_summary(path) -> dict:
    """What the split record keeps about an Imaris project, read from the .ims with h5py: the spot and
    surface counts, whether the Spots were created with region growing on, the manual Quality threshold
    the spots were kept with, the image shape and Imaris's own voxel size. Missing pieces give None
    rather than an error, so a project without spots or without creation parameters is recorded as
    unknown."""
    import h5py

    out = {"n_spots": None, "n_surfaces": None, "region_growing": None, "spot_quality_threshold": None,
           "shape_xyz": None, "imaris_voxel_um": None}
    with h5py.File(str(path), "r") as f:
        if _SPOTS in f:
            out["n_spots"] = int(f[_SPOTS].shape[0])
        if _SURF in f:
            out["n_surfaces"] = int(f[_SURF].shape[0])
        cp = _text_dataset(f, _POINTS_CP)
        if cp:
            out["region_growing"] = _xml_flag(cp, "mEnableRegionGrowing")
            block = re.search(r"<mSpotFilter>(.*?)</mSpotFilter>", cp, re.S)
            if block and 'mName="Quality"' in block.group(1):
                t = re.search(r'mLowerThresholdManualValue="([^"]+)"', block.group(1))
                if t:
                    out["spot_quality_threshold"] = float(t.group(1))
        info = f.get("DataSetInfo/Image")
        if info is not None:
            try:
                shape = [int(float(_attr_text(info, k))) for k in ("X", "Y", "Z")]
                ext_min = [float(_attr_text(info, f"ExtMin{i}")) for i in range(3)]
                ext_max = [float(_attr_text(info, f"ExtMax{i}")) for i in range(3)]
            except (TypeError, ValueError):
                pass
            else:
                out["shape_xyz"] = shape
                out["imaris_voxel_um"] = [round((b - a) / max(n, 1), 6) for a, b, n in zip(ext_min, ext_max, shape)]
    return out


def nd2_summary(path) -> dict:
    """Axis sizes and voxel size from the .nd2 metadata (no pixel data is read)."""
    import nd2

    with nd2.ND2File(str(path)) as f:
        vs = f.voxel_size()
        return {"sizes": {k: int(v) for k, v in f.sizes.items()},
                "voxel_um": [round(float(vs.x), 9), round(float(vs.y), 9), round(float(vs.z), 9)]}


# ----------------------------------------------------------------------------- the split


def build_split(real_manifest, ims_root, hash_cache: dict | None = None, registered_on: str | None = None,
                read_metadata: bool = True) -> dict:
    """Assemble the split record from the n2rad51 entries of a real manifest (their .nd2 paths) and
    the Imaris projects under ``ims_root`` (<ims_root>/<noHS_ims2|HS_ims2>/<stem>.ims). Every file is
    hashed unless ``hash_cache`` (see load_hash_cache) has its digest at the same size. With
    ``read_metadata`` the .nd2 sizes and voxel size and the Imaris summary are added."""
    runs = json.loads(Path(real_manifest).read_text(encoding="utf-8"))
    by_tag = {r["tag"]: r for r in runs if r.get("tag") in PREREGISTERED_ROLES}
    missing = sorted(set(PREREGISTERED_ROLES) - set(by_tag))
    if missing:
        raise ValueError(f"{real_manifest}: no entry for the registered tag(s) {missing}")
    gonads = []
    for tag, role in PREREGISTERED_ROLES.items():
        parts = tag_parts(tag)
        nd2_path = Path(by_tag[tag]["nd2"])
        stem = nd2_path.stem
        ims_path = Path(ims_root) / IMS_SUBDIR[parts["condition"]] / f"{stem}.ims"
        entry = {"tag": tag, "role": role, "condition": parts["condition"], "sex": parts["sex"],
                 "replicate": parts["replicate"], "stem": stem,
                 "nd2": _record_with_cache(nd2_path, hash_cache),
                 "ims": _record_with_cache(ims_path, hash_cache)}
        if read_metadata:
            entry["nd2"].update(nd2_summary(nd2_path))
            entry["ims"].update(read_ims_summary(ims_path))
        else:
            entry["ims"]["region_growing"] = None
        notes = list(NOTES.get(tag, []))
        rg = entry["ims"].get("region_growing")
        if rg is True:
            notes.append(REGION_GROWING_NOTE)
        elif rg is None:
            notes.append(REGION_GROWING_UNKNOWN_NOTE)
        entry["notes"] = notes
        gonads.append(entry)
    split = {
        "name": "n2rad51",
        "registered_on": registered_on or REGISTERED_ON,
        "protocol": "docs/BENCHMARK_PROTOCOL.md",
        "description": ("The 12 N2 RAD-51 gonads of 2025-11-05 (DAPI, SYP, RAD-51; no heat shock and heat shock; "
                        "hermaphrodites and males) with their Imaris annotations, split into train, val and test "
                        "before any benchmark result was computed. Roles are fixed; re-running "
                        "scripts/make_rad51_split.py only refreshes the file facts."),
        "frame": {
            "rule": ("Imaris coordinates are voxel centres scaled by Imaris's rounded voxel size (0.108 um in xy; "
                     "the .nd2 says 0.108333). Convert an Imaris coordinate to a voxel index with Imaris's own "
                     "voxel size, floor((coord - ExtMin) / imaris_voxel), and place it in the pipeline frame as "
                     "index times the .nd2 voxel size. Matching in raw micrometres drifts by almost a "
                     "micrometre across the 2600-pixel field."),
        },
        "roles": {role: sorted(g["tag"] for g in gonads if g["role"] == role) for role in ROLES},
        "gonads": gonads,
    }
    validate_split(split)
    return split


def validate_split(split: dict) -> None:
    """Raise ValueError when the split record is malformed: tags not unique or not n2rad51 tags, a role
    outside train/val/test, a file record without path, positive size and a 64-hex sha256, condition or
    sex disagreeing with the tag, or a ``roles`` header that disagrees with the gonad entries."""
    gonads = split.get("gonads")
    if not isinstance(gonads, list) or not gonads:
        raise ValueError("split has no 'gonads' list")
    tags = [g.get("tag") for g in gonads]
    if len(set(tags)) != len(tags):
        raise ValueError("split tags are not unique")
    for g in gonads:
        parts = tag_parts(g["tag"])
        if g.get("role") not in ROLES:
            raise ValueError(f"{g['tag']}: role {g.get('role')!r} is not one of {ROLES}")
        if g.get("condition") != parts["condition"] or g.get("sex") != parts["sex"]:
            raise ValueError(f"{g['tag']}: condition/sex disagree with the tag")
        for key in ("nd2", "ims"):
            rec = g.get(key)
            if not isinstance(rec, dict) or not rec.get("path"):
                raise ValueError(f"{g['tag']}: {key} record has no path")
            if not isinstance(rec.get("size_bytes"), int) or rec["size_bytes"] <= 0:
                raise ValueError(f"{g['tag']}: {key} size_bytes is not a positive integer")
            if not isinstance(rec.get("sha256"), str) or not _SHA256_RE.match(rec["sha256"]):
                raise ValueError(f"{g['tag']}: {key} sha256 is not a 64-digit hex string")
        if "stem" in g and Path(g["nd2"]["path"]).stem != g["stem"]:
            raise ValueError(f"{g['tag']}: stem {g['stem']!r} is not the .nd2 file name")
    header = split.get("roles")
    if header is not None:
        expect = {role: sorted(g["tag"] for g in gonads if g["role"] == role) for role in ROLES}
        got = {role: sorted(header.get(role, [])) for role in ROLES}
        if got != expect:
            raise ValueError("the 'roles' header disagrees with the gonad entries")


def load_split(path) -> dict:
    """Read and validate a split file."""
    split = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_split(split)
    return split


def split_tags(split: dict, role: str | None = None) -> list[str]:
    """Tags of the split, optionally those with one role, in file order."""
    return [g["tag"] for g in split["gonads"] if role is None or g["role"] == role]


def gonad(split: dict, tag: str) -> dict:
    for g in split["gonads"]:
        if g["tag"] == tag:
            return g
    raise KeyError(tag)


def check_split_files(split: dict, hashes: bool = False) -> list[str]:
    """Problems with the registered files: missing, a different size, and with ``hashes`` a different
    sha256 (which reads every file; about a minute for the 24 files when they are in the disk cache).
    An empty list means the data on disk is the data that was registered."""
    problems = []
    for g in split["gonads"]:
        for key in ("nd2", "ims"):
            rec = g[key]
            p = Path(rec["path"])
            if not p.is_file():
                problems.append(f"{g['tag']}: {key} missing: {rec['path']}")
                continue
            size = p.stat().st_size
            if size != rec["size_bytes"]:
                problems.append(f"{g['tag']}: {key} size {size} differs from the registered {rec['size_bytes']}")
                continue
            if hashes and file_sha256(p) != rec["sha256"]:
                problems.append(f"{g['tag']}: {key} sha256 differs from the registered digest")
    return problems


# ----------------------------------------------------------------------------- the manifest


def _exists_record(path, expected_size: int | None = None) -> dict:
    p = Path(path)
    exists = p.is_file()
    rec = {"path": posix(p), "exists": exists, "size_bytes": p.stat().st_size if exists else None}
    if expected_size is not None:
        rec["size_matches"] = (rec["size_bytes"] == expected_size) if exists else None
    return rec


def build_manifest(split: dict, results_root, truth_root, name: str | None = None,
                   truth_files: dict | None = None, split_path=None) -> dict:
    """One entry per gonad of the split: role, condition, sex, stem, the .nd2 and .ims (exists, size
    matches the registered size), the pipeline outputs under <results_root>/<tag>/ and the truth files
    under <truth_root>, each with its path and an exists flag, plus ``outputs_complete`` (nuclei table,
    label image and spots table all present), ``truth_complete`` (spots, surfaces and masks present),
    ``spot_truth_usable`` (false when the spots were exported with region growing, or when that is
    unknown) and the split's notes. Nothing is computed from the data."""
    patterns = dict(TRUTH_FILES, **(truth_files or {}))
    results_root, truth_root = Path(results_root), Path(truth_root)
    truth_names = sorted(p.name for p in truth_root.iterdir()) if truth_root.is_dir() else []
    entries = []
    for g in split["gonads"]:
        stem, tag = g["stem"], g["tag"]
        rdir = results_root / tag
        outputs = {k: _exists_record(rdir / pat.format(stem=stem)) for k, pat in OUTPUT_FILES.items()}
        truth = {k: _exists_record(truth_root / pat.format(stem=stem)) for k, pat in patterns.items()}
        found = [n for n in truth_names if n.startswith(f"{stem}__imaris_")]
        rg = g["ims"].get("region_growing")
        entries.append({
            "tag": tag, "role": g["role"], "condition": g["condition"], "sex": g["sex"], "stem": stem,
            "nd2": _exists_record(g["nd2"]["path"], g["nd2"].get("size_bytes")),
            "ims": {**_exists_record(g["ims"]["path"], g["ims"].get("size_bytes")), "region_growing": rg},
            "results_dir": posix(rdir),
            "outputs": outputs,
            "outputs_complete": all(outputs[k]["exists"] for k in REQUIRED_OUTPUTS),
            "truth": {**truth, "found": found},
            "truth_complete": all(truth[k]["exists"] for k in REQUIRED_TRUTH if k in truth),
            "spot_truth_usable": rg is False,
            "notes": list(g.get("notes", [])),
        })
    manifest = {
        "name": name or split.get("name", "benchmark"),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "protocol": split.get("protocol", "docs/BENCHMARK_PROTOCOL.md"),
        "split_file": ({"path": posix(split_path), "sha256": file_sha256(split_path),
                        "registered_on": split.get("registered_on")} if split_path else None),
        "results_root": posix(results_root),
        "truth_root": posix(truth_root),
        "truth_files": patterns,
        "roles": {role: [e["tag"] for e in entries if e["role"] == role] for role in ROLES},
        "summary": {
            "n_gonads": len(entries),
            "outputs_complete": [e["tag"] for e in entries if e["outputs_complete"]],
            "truth_complete": [e["tag"] for e in entries if e["truth_complete"]],
            "spot_truth_usable": [e["tag"] for e in entries if e["spot_truth_usable"]],
        },
        "gonads": entries,
    }
    return manifest


def write_manifest(manifest: dict, out_root) -> Path:
    """Write benchmarks/<name>/manifest.json under ``out_root`` and return its path."""
    out = Path(out_root) / manifest["name"] / "manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    return out


def today() -> str:
    return date.today().isoformat()
