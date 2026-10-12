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

* The additions (protocol changelog entry 9). ``tests/golden/rad51_additions.json`` lists gonads of
  another imaging session added to TRAIN after the split was registered, in a record of their own, so
  the split file stays byte-identical (``REGISTERED_SPLIT_SHA256``) and no registered role moves.
  ``ADDITION_ROLES`` fixes their roles in code (train only; an addition is never val or test and never
  stands in for a val gonad in a leave-one-gonad-out fold); ``build_additions`` assembles the record
  (scripts/make_rad51_additions.py), ``load_additions`` reads and validates it, ``check_addition_files``
  confirms the files are the recorded ones, and ``build_manifest(..., additions=...)`` appends them to a
  manifest without changing the registered entries. ``SEALED_NAME_MARKERS`` names files that are never
  read for the additions: the 2025-11-05 heat-shock male exports (the sealed test male's among them) and the
  lab's aggregate tables that compile the per-nucleus rows; a 2025-11-05 file is read only when it is a
  registered train or val project, by its stem behind the role guard (``guard_truth_read``,
  ``REGISTERED_STEMS``; protocol changelog entry 10).
* The spot truth of an addition (protocol changelog entries 9 (d) and 10): ``addition_spot_truth_usable``
  applies the export's own checks, then the decision of a changelog entry (``load_spot_truth_decisions``,
  tests/golden/rad51_additions_decisions.json), else the section 4.4 minimum (``SECTION_4_4_MIN``, read from
  the table of scripts/check_truth_assignment.py by ``load_truth_check``).
* The recorded spot sources (protocol changelog entry 11). ``tests/golden/rad51_spot_sources.json`` names, for a
  registered TRAIN gonad, an earlier Imaris statistics export of its own project whose spots replace the .ims's
  saved Spots object as its spot truth (n2rad51_HS_male_000: the curated 7,628 spots of its 2025-12-31 export,
  matched by position to the 772,485 candidates of the registered .ims). The split file does not change and no
  role moves (rule 4). ``SPOT_SOURCE_EXPORTS`` pins each export by name, size and sha256, so the sealed-name rule
  admits exactly that file (``spot_source_tag_of``, ``guard_export_read``, ``guard_truth_read``);
  ``load_spot_sources`` reads and validates the record; ``build_manifest(..., spot_sources=...)`` points the
  gonad's truth at the export from that source and decides its spot truth by ``spot_source_truth_usable``.

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
TRUTH_FILES = {   # exactly what `germquant imaris-truth` writes (src/germquant/validate/imaris_ims.py TRUTH_SUFFIXES)
    "spots": "{stem}__imaris_spots.csv",
    "nuclei": "{stem}__imaris_nuclei.csv",
    "labels": "{stem}__imaris_labels.tif",
    "params": "{stem}__imaris_params.json",
}
REQUIRED_TRUTH = ("spots", "nuclei", "labels")

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


# ----------------------------------------------------------------------------- the additions (entry 9)

ADDITIONS_NAME = "n2rad51_additions"
ADDITIONS_REGISTERED_ON = "2026-10-11"
ADDITIONS_CHANGELOG_ENTRY = 9
#: the session of the registered split (every gonad of PREREGISTERED_ROLES was imaged that day)
REGISTERED_SESSION = "2025-11-05"
#: sha256 of tests/golden/rad51_split.json as registered (and as it was when entry 9 added the additions): the
#: workstation checkout (CRLF line endings, the digest every manifest records) and the same bytes with LF line
#: endings (the git blob, a Linux checkout). The file changes only by the re-export of protocol section 4.1, with a
#: changelog line and new pinned digests here, in the same commit (protocol changelog entry 10 (g)(8));
#: tests/test_rad51_additions.py pins both, so any edit has to be deliberate and logged
REGISTERED_SPLIT_SHA256 = "955eeff52b4e953ca34520b2eae6a5a019520203e79967082a78eac4c991e25f"
REGISTERED_SPLIT_SHA256_LF = "989bfa08e4acd911004e433469979a2cfd52891aee240ac2f87db9183f0c9d85"
#: the additions of protocol changelog entry 9 and their roles: train only, fixed here so that no manifest or
#: record can move them (detector.gonad_role refuses a conflict); an addition leaves the benchmark only for a
#: defect on the truth side, recorded in the changelog
ADDITION_ROLES = {
    "n2rad51_20260105a_HS_male_001": "train",
    "n2rad51_20260105a_HS_male_002": "train",
    "n2rad51_20260105a_HS_male_003": "train",
    "n2rad51_20260105b_HS_male_001": "train",
}
ADDITION_TAG_RE = re.compile(
    r"^n2rad51_(?P<session>\d{8}[a-z])_(?P<condition>noHS|HS)_(?P<sex>herm|male)_(?P<replicate>\d{3})$")
#: lower-case markers of file names that are never read for the additions (protocol changelog entries 9 and 10):
#: no test-gonad file, no 2025-11-05 heat-shock male export and no aggregate of the per-nucleus tables (the
#: decision-1 export check of entry 9 lists them). A file of the 2025-11-05 session is refused unless it is a
#: registered project read by its registered stem behind the role guard (``guard_truth_read``): the eight train
#: and val projects were read that way for entry 9 (d)(3) and (d)(4); or the one recorded earlier export of a
#: registered train gonad that protocol changelog entry 11 admits by its pinned name and size
#: (``SPOT_SOURCE_EXPORTS``), behind the role guard
SEALED_NAME_MARKERS = ("20251105", "temp_compiledtdata", "compiled_data", "master_data", "rad51_stats", ".prism")
#: the stems of the registered split (tests/golden/rad51_split.json; tests/test_rad51_additions.py checks them
#: against the file), so a truth-side read can tell a registered project, which the role guard admits or refuses,
#: from any other 2025-11-05 file, which is sealed
REGISTERED_STEMS = {
    "n2rad51_noHS_herm_000": "20251105_N2_nohs_HERM _",
    "n2rad51_noHS_herm_001": "20251105_N2_nohs_HERM _001",
    "n2rad51_noHS_herm_002": "20251105_N2_nohs_HERM _002",
    "n2rad51_noHS_male_000": "20251105_N2_nohs_MALE _",
    "n2rad51_noHS_male_001": "20251105_N2_nohs_MALE _001",
    "n2rad51_noHS_male_002": "20251105_N2_nohs_MALE _002",
    "n2rad51_HS_herm_000": "20251105_n2_hs_HERM _",
    "n2rad51_HS_herm_001": "20251105_n2_hs_HERM _001",
    "n2rad51_HS_herm_002": "20251105_n2_hs_HERM _002",
    "n2rad51_HS_male_000": "20251105_n2_hs_MALE _",
    "n2rad51_HS_male_001": "20251105_n2_hs_MALE _001",
    "n2rad51_HS_male_002": "20251105_n2_hs_MALE _002",
}
SPOT_TRUTH_SOURCES = ("ims", "xlsx")
SPOT_TRUTH_STATUS = ("usable", "pending")
ADDITION_FILE_KEYS = ("nd2", "ims", "xlsx")

#: protocol changelog entry 11: the registered train gonads whose spot truth comes from a recorded earlier Imaris
#: statistics export of their own project instead of the .ims's saved Spots object (the record is
#: tests/golden/rad51_spot_sources.json). Each export is pinned here by file name, size and sha256: the guard admits
#: exactly that file (by name and size, before anything is read; the sha256 is checked when it is read), behind the
#: role guard. n2rad51_HS_male_000's registered .ims holds 772,485 unthresholded, unclassified candidates of constant
#: radius; its 2025-12-31 export (MALE_3 of the three, identified by name, size and date) holds the curated 7,628,
#: region growing off by the section 4.1 rule. The test male's exports (MALE_1) carry another name and size and stay
#: sealed; HS_male_002's (MALE_2) is region-grown and is not a spot source.
SPOT_SOURCES_CHANGELOG_ENTRY = 11
SPOT_SOURCE_KIND = "earlier_export_subset"
SPOT_SOURCE_EXPORTS = {
    "n2rad51_HS_male_000": {"name": "20251105_n2_hs_MALE_3.xlsx", "size_bytes": 15512227,
                            "sha256": "0fbc41b7dc5dbf0a1ef0c87f367451635c0f7dfea1302e6dfc53b614082aca52"},
}
SPOT_SOURCE_NOTE = (
    "Spot truth from a recorded earlier export (protocol changelog entry 11, tests/golden/rad51_spot_sources.json): "
    "the curated spots of the project's 2025-12-31 statistics export, matched by position to the spots of the "
    "registered .ims, Class A and B by the project's own classification rule. The .ims, its surfaces and the nucleus "
    "truth are the registered ones.")

# the protocol's text (sections 4.1 and 4.2) as values, for the creation-parameter comparison of entry 9 (d)(4)
PROTOCOL_TEXT_PARAMETERS = {
    "spots.source_channel_index": 2, "spots.estimated_diameter_um": [0.4, 0.4, 0.9],
    "spots.background_subtraction": True, "spots.roi_enabled": False, "spots.region_growing": False,
    "surfaces.source_channel_index": 0, "surfaces.grain_um": 0.4, "surfaces.background_subtraction": True,
    "surfaces.largest_sphere_um": 2.5, "surfaces.automatic_threshold": False, "surfaces.seed_split": True,
    "surfaces.seed_diameter_um": 2.5, "surfaces.voxel_filter_enabled": True, "surfaces.roi_enabled": True,
}
# numbers that vary per gonad by design; compared with the range of the registered projects, not with one value
PER_GONAD_RANGE_KEYS = ("spots.quality_threshold", "spots.class_threshold_um", "surfaces.manual_threshold",
                        "surfaces.voxel_filter_min", "surfaces.voxel_filter_max", "surfaces.class_threshold")


def canonical_json(obj) -> str:
    """Key-sorted compact JSON (the form a sha256 of a mapping is taken of)."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def roles_sha256(roles: dict | None = None) -> str:
    """sha256 of the canonical JSON of the registered roles (``PREREGISTERED_ROLES`` by default). The additions
    record carries it, so a record cannot be read against a split whose roles moved."""
    return hashlib.sha256(canonical_json(PREREGISTERED_ROLES if roles is None else roles).encode("utf-8")).hexdigest()


def split_file_unchanged(path) -> bool:
    """True when ``path`` holds the registered split file byte for byte (with CRLF or with LF line endings)."""
    b = Path(path).read_bytes()
    return (hashlib.sha256(b).hexdigest() == REGISTERED_SPLIT_SHA256
            or lf_sha256(b) == REGISTERED_SPLIT_SHA256_LF)


def lf_sha256(data: bytes) -> str:
    """sha256 of ``data`` with CRLF line endings turned into LF (the git blob of a text file)."""
    return hashlib.sha256(data.replace(b"\r\n", b"\n")).hexdigest()


def addition_tag_parts(tag: str) -> dict:
    """session (yyyymmdd plus a set letter), condition, sex and replicate of an addition's tag."""
    m = ADDITION_TAG_RE.match(str(tag))
    if not m:
        raise ValueError(f"{tag!r} is not an addition tag (n2rad51_<yyyymmdd><set letter>_<noHS|HS>_<herm|male>_<nnn>)")
    return m.groupdict()


def addition_role(tag: str) -> str | None:
    """The role entry 9 gives an addition (train), None for any other tag."""
    return ADDITION_ROLES.get(str(tag))


def is_sealed_path(path) -> bool:
    """True when a file name or folder of ``path`` carries a sealed marker (``SEALED_NAME_MARKERS``)."""
    text = str(path).replace("\\", "/").lower()
    return any(m in text for m in SEALED_NAME_MARKERS)


def refuse_sealed(path) -> None:
    """Raise PermissionError for a sealed path before anything opens it."""
    if is_sealed_path(path):
        raise PermissionError(f"refused: {path} is sealed (protocol changelog entries 9 and 10: no test-gonad file, "
                              "no 2025-11-05 heat-shock male export and no aggregate of the per-nucleus tables is "
                              "read; a registered train or val project is read only by its registered stem, behind "
                              "the role guard)")


def registered_tag_of(path) -> str | None:
    """The tag of the registered gonad whose stem is the file name of ``path`` without its extension (letter case
    and surrounding spaces ignored), None for any other file."""
    stem = Path(str(path).replace("\\", "/")).stem.strip().lower()
    for tag, s in REGISTERED_STEMS.items():
        if s.strip().lower() == stem:
            return tag
    return None


def spot_source_tag_of(path) -> str | None:
    """The registered gonad whose pinned earlier export (``SPOT_SOURCE_EXPORTS``, protocol changelog entry 11) is
    ``path``: the file name equals the pinned name (letter case ignored) and the file on disk has the pinned size.
    Only the name and the file size are looked at, so nothing is read; None for any other file (a missing file, a
    file of another size, any other name)."""
    p = Path(str(path).replace("\\", "/"))
    for tag, rec in SPOT_SOURCE_EXPORTS.items():
        if p.name.lower() == rec["name"].lower():
            try:
                return tag if p.is_file() and p.stat().st_size == rec["size_bytes"] else None
            except OSError:
                return None
    return None


def guard_export_read(path) -> str | None:
    """Before an Imaris statistics export is opened (``imaris_xlsx``): the pinned earlier export of a registered gonad
    (``spot_source_tag_of``) passes the role guard (it is a train gonad's, or the guard raises); any other file must
    not be sealed (``refuse_sealed``). Returns the registered tag of a pinned export, None otherwise."""
    tag = spot_source_tag_of(path)
    if tag is None:
        refuse_sealed(path)
        return None
    from germquant.benchmark.detector import guard_test

    guard_test([tag])
    return tag


def guard_truth_read(paths, *, allow_test: bool = False) -> list[str]:
    """Before a truth-side read (the frame check, the xlsx spot source or the recorded earlier export of ``germquant
    imaris-truth``): a file of a registered gonad passes the role guard (train and val; test only with
    ``allow_test``, the pre-flight of docs/FINAL_SCORING.md), the pinned earlier export of a registered gonad
    (``spot_source_tag_of``, protocol changelog entry 11) counts as that gonad's file, and any other file must not
    be sealed (``refuse_sealed``). None entries are skipped. Returns the registered tags found. Raises before
    anything is opened."""
    from germquant.benchmark.detector import guard_test

    tags = []
    for p in paths:
        if p is None:
            continue
        tag = registered_tag_of(p) or spot_source_tag_of(p)
        if tag is None:
            refuse_sealed(p)
        elif tag not in tags:
            tags.append(tag)
    guard_test(tags, allow_test=allow_test)
    return tags


# ---- the Imaris and nd2 facts of a project ----

def _xml_attr(text: str | None, name: str) -> str | None:
    m = re.search(rf'\b{name}="([^"]*)"', text or "")
    return None if m is None else m.group(1)


def _xml_bool(text: str | None, name: str) -> bool | None:
    v = _xml_attr(text, name)
    return None if v is None else v.lower() == "true"


def _xml_float(text: str | None, name: str) -> float | None:
    v = _xml_attr(text, name)
    try:
        return None if v is None else float(v)
    except ValueError:
        return None


def _filters(text: str | None, block: str) -> list[dict]:
    """The statistics filters inside ``<block>...</block>`` of a creation-parameter text: statistic name,
    lower and upper threshold (enabled flag, manual flag, value)."""
    m = re.search(rf"<{block}>(.*?)</{block}>", text or "", re.S)
    out = []
    if not m:
        return out
    for f in re.finditer(r"<bpStatisticsFilter\s+([^>]*)>(.*?)</bpStatisticsFilter>", m.group(1), re.S):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', f.group(1)))
        name = re.search(r'mName="([^"]*)"', f.group(2))

        def num(k):
            try:
                return float(attrs[k])
            except (KeyError, ValueError):
                return None

        out.append({"statistic": name.group(1) if name else None,
                    "lower_enabled": attrs.get("mLowerThresholdEnable", "").lower() == "true",
                    "lower_manual": attrs.get("mLowerThresholdManual", "").lower() == "true",
                    "lower": num("mLowerThresholdManualValue"),
                    "upper_enabled": attrs.get("mUpperThresholdEnable", "").lower() == "true",
                    "upper": num("mUpperThresholdManualValue")})
    return out


def _classification(text: str | None) -> dict:
    """The one-dimensional class split of a creation-parameter text: statistic and the Class A upper bound."""
    m = re.search(r"<bpCObjectsClassifyFilter1D>(.*?)</bpCObjectsClassifyFilter1D>", text or "", re.S)
    if not m:
        return {"statistic": None, "threshold": None}
    stat = re.search(r'mName="([^"]*)"', m.group(1))
    thr = re.search(r'<cThresholdValue mClassName="0"[^>]*mThresholdHigh="([^"]*)"', m.group(1))
    return {"statistic": stat.group(1) if stat else None, "threshold": float(thr.group(1)) if thr else None}


def _roi_boxes(text: str | None) -> list[dict]:
    out = []
    for m in re.finditer(r"<bpRegionOfInterest\s+([^>]*?)/?>", text or ""):
        a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        out.append({"min_xyz_voxel": [int(float(v)) for v in a.get("mMinXYZ", "").split()],
                    "max_xyz_voxel": [int(float(v)) for v in a.get("mMaxXYZ", "").split()]})
    return out


def creation_parameter_summary(points_cp: str | None, surfaces_cp: str | None) -> dict:
    """The creation parameters entry 9 (d)(4) compares, from the Spots and the Surfaces creation-parameter
    texts of an Imaris project (``{"spots": {...}, "surfaces": {...}}``; None where a text lacks a value)."""
    spot_filters = _filters(points_cp, "mSpotFilter")
    quality = [f for f in spot_filters if f["statistic"] == "Quality"]
    diam = _xml_attr(points_cp, "mEstimatedDiameter")
    spots = {
        "source_channel_index": int(_xml_float(points_cp, "mSourceChannelIndex")) if _xml_attr(points_cp, "mSourceChannelIndex") else None,
        "estimated_diameter_um": [round(float(v), 6) for v in diam.split()] if diam else None,
        "background_subtraction": _xml_bool(points_cp, "mBackgroundSubtraction"),
        "roi_enabled": _xml_bool(points_cp, "mEnableRegionOfInterest"),
        "region_growing": _xml_bool(points_cp, "mEnableRegionGrowing"),
        "quality_threshold": quality[0]["lower"] if quality else None,
        "spot_filters": [f["statistic"] for f in spot_filters],
        "class_statistic": _classification(points_cp)["statistic"],
        "class_threshold_um": _classification(points_cp)["threshold"],
        "roi_boxes": _roi_boxes(points_cp),
    }
    voxel_filter = [f for f in _filters(surfaces_cp, "mSurfacesFilter") if f["statistic"] == "Number of Voxels"]
    vf = voxel_filter[0] if voxel_filter else None
    surfaces = {
        "source_channel_index": int(_xml_float(surfaces_cp, "mSourceChannelIndex")) if _xml_attr(surfaces_cp, "mSourceChannelIndex") else None,
        "smoothing": _xml_bool(surfaces_cp, "mEnableSmooth"),
        "grain_um": _xml_float(surfaces_cp, "mSurfaceGrainSize"),
        "background_subtraction": _xml_bool(surfaces_cp, "mEnableEliminateBackground"),
        "largest_sphere_um": _xml_float(surfaces_cp, "mDiameterOfLargestSphere"),
        "automatic_threshold": _xml_bool(surfaces_cp, "mEnableAutomaticThreshold"),
        "manual_threshold": _xml_float(surfaces_cp, "mManualThresholdValue"),
        "seed_split": _xml_bool(surfaces_cp, "mEnableRegionGrowing"),
        "seed_diameter_um": _xml_float(surfaces_cp, "mRegionGrowingEstimatedDiameter"),
        "voxel_filter_enabled": bool(vf and (vf["lower_enabled"] or vf["upper_enabled"])),
        "voxel_filter_min": vf["lower"] if vf and vf["lower_enabled"] else None,
        "voxel_filter_max": vf["upper"] if vf and vf["upper_enabled"] else None,
        "surface_filters": [f["statistic"] for f in _filters(surfaces_cp, "mSurfacesFilter")],
        "roi_enabled": _xml_bool(surfaces_cp, "mEnableRegionOfInterest"),
        "roi_boxes": _roi_boxes(surfaces_cp),
        "class_statistic": _classification(surfaces_cp)["statistic"],
        "class_threshold": _classification(surfaces_cp)["threshold"],
    }
    return {"spots": spots, "surfaces": surfaces}


def _flat(summary: dict) -> dict:
    return {f"{obj}.{k}": v for obj, d in summary.items() for k, v in d.items() if k not in ("roi_boxes",)}


def compare_creation_parameters(summary: dict, registered: dict[str, dict]) -> list[dict]:
    """Where a project's creation parameters (``creation_parameter_summary``) differ from the protocol's text
    (sections 4.1 and 4.2, ``PROTOCOL_TEXT_PARAMETERS``) and from the registered projects (``registered``: tag
    -> summary; entry 9 uses the eight train and val projects). Per-gonad numbers (``PER_GONAD_RANGE_KEYS``)
    are compared with the registered range; every other key with the set of registered values. A difference is
    recorded, not a defect."""
    mine = _flat(summary)
    reg = {t: _flat(s) for t, s in registered.items()}
    out = []
    for key, value in mine.items():
        rec = {"parameter": key, "value": value}
        if key in PROTOCOL_TEXT_PARAMETERS and value != PROTOCOL_TEXT_PARAMETERS[key]:
            rec["protocol_text"] = PROTOCOL_TEXT_PARAMETERS[key]
        vals = [r.get(key) for r in reg.values()]
        if key in PER_GONAD_RANGE_KEYS:
            nums = [v for v in vals if isinstance(v, (int, float)) and not isinstance(v, bool)]
            if value is None:
                if reg and not any(v is None for v in vals):
                    rec["registered_range"] = [min(nums), max(nums)] if nums else None
            elif nums:
                lo, hi = min(nums), max(nums)
                if isinstance(value, bool) or not lo <= value <= hi:
                    rec["registered_range"] = [lo, hi]
        else:
            seen = sorted({json.dumps(v, sort_keys=True) for v in vals})
            if reg and json.dumps(value, sort_keys=True) not in seen:
                rec["registered_values"] = [json.loads(s) for s in seen]
        if len(rec) > 2:
            out.append(rec)
    return out


def read_ims_facts(path) -> dict:
    """What the additions record keeps about an Imaris project (entry 9 (b)): the summary of ``read_ims_summary``
    plus the raw extent text (six significant digits; ``frame`` in the record explains why it matters), the
    per-axis Imaris voxel, channel names, Imaris version, the spot source channel, the stored spot statistics,
    whether a Masked channel is present and the creation parameters (``creation_parameter_summary``)."""
    import h5py

    refuse_sealed(path)
    out = read_ims_summary(path)
    with h5py.File(str(path), "r") as f:
        info = f["DataSetInfo/Image"]
        out["ext_min_xyz_text"] = [_attr_text(info, f"ExtMin{i}") for i in range(3)]
        out["ext_max_xyz_text"] = [_attr_text(info, f"ExtMax{i}") for i in range(3)]
        lo = [float(v) for v in out["ext_min_xyz_text"]]
        hi = [float(v) for v in out["ext_max_xyz_text"]]
        out["imaris_voxel_um"] = [round((b - a) / max(n, 1), 7) for a, b, n in zip(lo, hi, out["shape_xyz"])]
        names = []
        for k in sorted((k for k in f["DataSetInfo"].keys() if k.startswith("Channel ")), key=lambda k: int(k.split()[-1])):
            names.append(_attr_text(f["DataSetInfo"][k], "Name"))
        out["channel_names"] = names
        out["masked_channel_present"] = any("mask" in (n or "").lower() for n in names)
        out["imaris_version"] = _attr_text(f["DataSetInfo/Imaris"], "Version") if "DataSetInfo/Imaris" in f else None
        pcp = _text_dataset(f, _POINTS_CP)
        scp = _text_dataset(f, "Scene8/Content/MegaSurfaces0/CreationParameters")
        summ = creation_parameter_summary(pcp, scp)
        idx = summ["spots"]["source_channel_index"]
        out["spot_source_channel"] = {"index": idx, "name": names[idx] if idx is not None and idx < len(names) else None}
        stats = set()
        g = f.get("Scene8/Content/Points0")
        if g is not None and "StatisticsType" in g:
            cats = {int(r["ID"]): r["CategoryName"].decode("latin1") for r in g["Category"][:]} if "Category" in g else {}
            for r in g["StatisticsType"][:]:
                if cats.get(int(r["ID_Category"])) != "Overall":
                    stats.add(r["Name"].decode("latin1").rstrip("\x00"))
        out["spot_statistics"] = sorted(stats)
    out["creation_parameters"] = summ
    return out


def nd2_facts(path) -> dict:
    """The .nd2 header facts of the additions record: axis sizes, voxel size (x, y, z), channel names, dtype,
    the camera and each channel's exposure (from the Nikon description). No pixel is read."""
    import nd2

    from ..io.acquisition import parse_description, read_description

    refuse_sealed(path)
    out = nd2_summary(path)
    with nd2.ND2File(str(path)) as f:
        out["channel_names"] = [str(c.channel.name) for c in (f.metadata.channels or [])]
        out["dtype"] = str(f.dtype)
    try:
        desc = read_description(path)
        cam = re.search(r"Camera Name:\s*([^\r\n]+)", desc)
        out["camera"] = cam.group(1).strip() if cam else None
        per = parse_description(desc)
        out["exposure_ms"] = {name: rec.get("exposure_ms") for name, rec in per.items()}
    except Exception as e:  # noqa: BLE001 - the description is a convenience, not a fact the record needs
        out["camera"], out["exposure_ms"] = None, {"error": f"{type(e).__name__}: {e}"}
    return out


# ---- the record ----

def _copy_manifest_index(path) -> dict:
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    return {posix(r["dest"]): r for r in doc.get("files", [])}


def _addition_file(local, copies: dict, hash_files: bool) -> dict:
    """{path, nas_path, nas_mtime, size_bytes, sha256} of one copied file, its size and sha256 re-verified
    against the copy manifest (ValueError on a difference)."""
    refuse_sealed(local)
    p = Path(local)
    c = copies.get(posix(p))
    if c is None:
        raise ValueError(f"{p}: not in the copy manifest")
    refuse_sealed(c.get("source", ""))
    size = p.stat().st_size
    if size != int(c["size"]):
        raise ValueError(f"{p}: {size} bytes, the copy manifest says {c['size']}")
    sha = file_sha256(p) if hash_files else c["sha256"]
    if sha != c["sha256"]:
        raise ValueError(f"{p}: sha256 {sha} is not the copy manifest's {c['sha256']}")
    return {"path": posix(p), "nas_path": c.get("source"), "nas_mtime": c.get("source_mtime"), "size_bytes": size,
            "sha256": sha, "sha256_verified_against_copy_manifest": bool(hash_files)}


def build_additions(specs: list[dict], copy_manifest, *, registered_split_path, decided_by: dict,
                    registered_params: dict | None = None, not_added: list | None = None,
                    hash_files: bool = True, read_metadata: bool = True,
                    registered_on: str = ADDITIONS_REGISTERED_ON) -> dict:
    """The additions record (entry 9 (b)). ``specs``: one dict per gonad with ``tag``, ``set`` (the NAS folder),
    ``genotype``, ``stem``, local ``nd2``, ``ims`` and optional ``xlsx`` paths, ``xlsx_use``, ``spot_truth``
    ({source, status}), ``nucleus_truth`` and ``notes``. Every file is hashed and compared with the copy manifest;
    with ``read_metadata`` the nd2 header and the Imaris facts are read and the creation parameters compared with
    ``registered_params`` (tag -> ``creation_parameter_summary`` of the registered train and val projects)."""
    copies = _copy_manifest_index(copy_manifest)
    split_sha = file_sha256(registered_split_path)
    gonads = []
    for s in specs:
        tag = s["tag"]
        parts = addition_tag_parts(tag)
        g = {"tag": tag, "role": ADDITION_ROLES.get(tag), "session": s["session"], "set": s["set"],
             "genotype": s.get("genotype", "N2"), "condition": parts["condition"], "sex": parts["sex"],
             "replicate": parts["replicate"], "stem": s["stem"]}
        for key in ADDITION_FILE_KEYS:
            if s.get(key):
                g[key] = _addition_file(s[key], copies, hash_files)
        if read_metadata:
            g["nd2"].update(nd2_facts(g["nd2"]["path"]))
            g["ims"].update(read_ims_facts(g["ims"]["path"]))
        if "xlsx" in g:
            g["xlsx"]["use"] = s.get("xlsx_use")
            for k in ("overall_n_spots", "overall_n_surfaces", "overall_source"):
                if k in s:
                    g["xlsx"][k] = s[k]
        g["spot_truth"] = dict(s["spot_truth"])
        g["nucleus_truth"] = s.get("nucleus_truth")
        if read_metadata and registered_params:
            g["creation_parameter_differences"] = compare_creation_parameters(g["ims"]["creation_parameters"],
                                                                             registered_params)
        g["notes"] = list(s.get("notes", []))
        gonads.append(g)
    rec = {
        "name": ADDITIONS_NAME,
        "registered_on": registered_on,
        "protocol": "docs/BENCHMARK_PROTOCOL.md",
        "changelog_entry": ADDITIONS_CHANGELOG_ENTRY,
        "description": ("Gonads added to TRAIN after the split of tests/golden/rad51_split.json was registered "
                        "(protocol changelog entry 9). The split file and its roles are unchanged; these gonads are "
                        "never val or test, never stand in for a val gonad, and leave the benchmark only for a "
                        "defect on the truth side, recorded in the changelog."),
        "decided_by": decided_by,
        "registered_split": {"path": "tests/golden/rad51_split.json", "sha256": split_sha,
                             "sha256_lf": lf_sha256(Path(registered_split_path).read_bytes()),
                             "roles_sha256": roles_sha256(),
                             "unchanged": split_file_unchanged(registered_split_path)},
        "rules": [
            "train only: an addition is never val or test (ADDITION_ROLES in germquant.validate.benchmark_split)",
            "an addition never stands in for a held-out val gonad in a leave-one-gonad-out fold",
            "an addition leaves the benchmark only for a defect on the truth side, recorded in the changelog",
            "nothing runs on these files before this record and changelog entry 9 are committed",
            "no 2025-11-05 file and no aggregate of the 2025-11-05 per-nucleus tables is read (SEALED_NAME_MARKERS)",
        ],
        "frame": {
            "rule": ("As in the split record: floor((coordinate - ExtMin) / imaris_voxel) per axis, placed at index "
                     "times the .nd2 voxel."),
            "per_axis_voxel_note": ("Imaris stores the extents as text with six significant digits. On these "
                                    "2720-pixel fields at stage x above 10000 um the x extent reads 293.8 um "
                                    "(0.10801 um per voxel) where y reads 293.76 um (0.108 um), and ExtMin x is known "
                                    "only to +-0.05 um. The frame used for each gonad is decided by the test of entry "
                                    "9 (d)(3) and recorded in its truth params JSON (frame_check)."),
        },
        "copy_manifest": {"path": posix(copy_manifest), "sha256": file_sha256(copy_manifest)},
        "roles": {r: sorted(g["tag"] for g in gonads if g["role"] == r) for r in ROLES},
        "gonads": gonads,
        "not_added": list(not_added or []),
    }
    validate_additions(rec)
    return rec


def validate_additions(rec: dict) -> None:
    """Raise ValueError when the additions record is malformed or would move a registered role: a tag that is
    not an addition tag, a registered tag, a tag outside ``ADDITION_ROLES`` or with another role than it gives
    (so never val or test), a duplicate, a file record without path, positive size and 64-hex sha256, a stem
    that is not the .nd2 name, condition or sex disagreeing with the tag, a sealed path, a spot-truth source or
    status outside the allowed ones, a roles header that disagrees with the entries, or a record made against
    other registered roles (``registered_split.roles_sha256``)."""
    gonads = rec.get("gonads")
    if not isinstance(gonads, list) or not gonads:
        raise ValueError("the additions record has no 'gonads' list")
    reg = rec.get("registered_split") or {}
    if reg.get("roles_sha256") != roles_sha256():
        raise ValueError("the additions record was made against other registered roles (roles_sha256 differs)")
    tags = [g.get("tag") for g in gonads]
    if len(set(tags)) != len(tags):
        raise ValueError("addition tags are not unique")
    for g in gonads:
        tag = g.get("tag")
        if tag in PREREGISTERED_ROLES:
            raise ValueError(f"{tag}: a registered gonad cannot be an addition (rule 4: the split never changes)")
        parts = addition_tag_parts(tag)
        if tag not in ADDITION_ROLES:
            raise ValueError(f"{tag}: not an addition of entry 9 (ADDITION_ROLES)")
        if g.get("role") != ADDITION_ROLES[tag]:
            raise ValueError(f"{tag}: role {g.get('role')!r}; entry 9 makes it {ADDITION_ROLES[tag]!r} "
                             "(an addition is never val or test)")
        if g.get("condition") != parts["condition"] or g.get("sex") != parts["sex"]:
            raise ValueError(f"{tag}: condition/sex disagree with the tag")
        for key in ADDITION_FILE_KEYS:
            r = g.get(key)
            if r is None and key == "xlsx":
                continue
            if not isinstance(r, dict) or not r.get("path"):
                raise ValueError(f"{tag}: {key} record has no path")
            if not isinstance(r.get("size_bytes"), int) or r["size_bytes"] <= 0:
                raise ValueError(f"{tag}: {key} size_bytes is not a positive integer")
            if not isinstance(r.get("sha256"), str) or not _SHA256_RE.match(r["sha256"]):
                raise ValueError(f"{tag}: {key} sha256 is not a 64-digit hex string")
            for k in ("path", "nas_path"):
                if r.get(k) and is_sealed_path(r[k]):
                    raise ValueError(f"{tag}: {key} {k} {r[k]} is sealed")
        if Path(g["nd2"]["path"]).stem != g.get("stem"):
            raise ValueError(f"{tag}: stem {g.get('stem')!r} is not the .nd2 file name")
        st = g.get("spot_truth") or {}
        if st.get("source") not in SPOT_TRUTH_SOURCES or st.get("status") not in SPOT_TRUTH_STATUS:
            raise ValueError(f"{tag}: spot_truth needs a source in {SPOT_TRUTH_SOURCES} and a status in "
                             f"{SPOT_TRUTH_STATUS}")
        if st["source"] == "xlsx" and not g.get("xlsx"):
            raise ValueError(f"{tag}: spot truth from the xlsx needs an xlsx record")
    header = rec.get("roles")
    if header is not None:
        expect = {role: sorted(g["tag"] for g in gonads if g["role"] == role) for role in ROLES}
        got = {role: sorted(header.get(role, [])) for role in ROLES}
        if got != expect:
            raise ValueError("the 'roles' header disagrees with the gonad entries")


def load_additions(path) -> dict:
    """Read and validate an additions record."""
    rec = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_additions(rec)
    return rec


def check_addition_files(rec: dict, hashes: bool = False) -> list[str]:
    """Problems with the recorded files of the additions (missing, another size, with ``hashes`` another
    sha256); an empty list means the files on disk are the recorded ones."""
    problems = []
    for g in rec["gonads"]:
        for key in ADDITION_FILE_KEYS:
            r = g.get(key)
            if r is None:
                continue
            p = Path(r["path"])
            if not p.is_file():
                problems.append(f"{g['tag']}: {key} missing: {r['path']}")
                continue
            if p.stat().st_size != r["size_bytes"]:
                problems.append(f"{g['tag']}: {key} size {p.stat().st_size} differs from the recorded {r['size_bytes']}")
                continue
            if hashes and file_sha256(p) != r["sha256"]:
                problems.append(f"{g['tag']}: {key} sha256 differs from the recorded digest")
    return problems


#: the section 4.4 minimum of entry 9 (d)(5) and (e) V7, read as a minimum by entry 10 (b): the lowest agreement of
#: the registered train and val gonads (entry 4 (a)); the same number as imaris_ims.TRUTH_ASSIGNMENT_MIN
SECTION_4_4_MIN = 0.9935
#: the xlsx criterion that is the section 4.4 check (entry 9 (e)); the others are checks of the export itself
XLSX_SECTION_4_4_CRITERION = "V7_truth_assignment"


def addition_spot_truth_usable(g: dict, truth_params: dict | None, truth_check: dict | None = None,
                               decision: dict | None = None) -> tuple[bool, str]:
    """Whether an addition's spot truth may be scored, and why. The rules, in order (protocol changelog entries 9
    (d), (e) and 10 (c)):

    1. no truth params (no truth export yet): not usable;
    2. the export itself: from the .ims, the record calls the spot truth usable and region growing is off; from
       the xlsx (entry 9 (e)), every criterion of its verification but the section 4.4 one (V7) passes; and in
       both cases the export ran the frame check of entry 9 (d)(3) and it chose a frame (an export made without
       ``--frame-check`` is not usable). A decision cannot override this;
    3. ``decision`` (an entry of the decisions file of a later changelog entry, ``load_spot_truth_decisions``):
       its ``usable`` and its reason;
    4. otherwise section 4.4: ``truth_check`` (this gonad's row of scripts/check_truth_assignment.py, at least
       ``agreement``; for an xlsx source V7's own agreement when no row is given) must reach ``SECTION_4_4_MIN``.
       A failed or missing check leaves the spot truth unusable until a changelog entry decides it."""
    if truth_params is None:
        return False, "no truth export yet"
    st = g.get("spot_truth") or {}
    params = truth_params
    agreement = None if truth_check is None else truth_check.get("agreement")
    if st.get("source") == "xlsx":
        ver = ((params.get("spots_source") or {}).get("verification") or {})
        crit = {k: v for k, v in ver.items() if k.startswith("V")}
        if not crit:
            return False, "xlsx spots not verified (entry 9 (e))"
        failed = sorted(k for k, v in crit.items() if not (v or {}).get("passes"))
        if any(k != XLSX_SECTION_4_4_CRITERION for k in failed):
            return False, f"xlsx spots failed {', '.join(k.split('_')[0] for k in failed)} of entry 9 (e)"
        if agreement is None:
            agreement = (crit.get(XLSX_SECTION_4_4_CRITERION) or {}).get("agreement")
    else:
        if g.get("ims", {}).get("region_growing") is not False:
            return False, "spot truth unusable (Imaris region growing, or unknown)"
        if st.get("status") != "usable":
            return False, f"spot truth status {st.get('status')!r}"
    fc = params.get("frame_check")
    if fc is None or fc.get("chosen") is None:
        return False, "frame check of entry 9 (d)(3) missing or failed"
    return _decision_or_section_4_4(agreement, decision, "entry 9 (d)(5)")


def _decision_or_section_4_4(agreement, decision: dict | None, rule: str) -> tuple[bool, str]:
    """The last two rules of the spot-truth decision: a changelog decision, else a section 4.4 agreement of at least
    ``SECTION_4_4_MIN`` (``rule`` names where the criterion was declared)."""
    if decision is not None:
        return bool(decision["usable"]), f"changelog entry {decision['entry']}: {decision['reason']}"
    try:
        a = float(agreement)
    except (TypeError, ValueError):
        a = float("nan")
    if not a >= SECTION_4_4_MIN:
        return False, (f"section 4.4 below 99.35 percent or not run ({rule}); needs a changelog entry"
                       + ("" if agreement is None else f" (agreement {a:.4f})"))
    return True, f"section 4.4 agreement {a:.4f}, at least {SECTION_4_4_MIN} ({rule})"


# ---- the recorded earlier exports of registered train gonads (protocol changelog entry 11) ----

#: the criterion of the export-subset verification (``imaris_ims``, entry 11) that is the section 4.4 check; the
#: others (E1 to E9) are checks of the export itself
EXPORT_SECTION_4_4_CRITERION = "E10_truth_assignment"
SPOT_SOURCE_RULE = "entry 11"


def validate_spot_sources(rec: dict, split: dict | None = None) -> None:
    """Raise ValueError when the spot-sources record (tests/golden/rad51_spot_sources.json, protocol changelog entry
    11) is malformed or names anything but the pinned export of a registered train gonad: the record of a changelog
    entry 11 or later, made against the registered roles; every source a registered TRAIN gonad (never val or test;
    rule 4: no role changes) with its registered stem, the kind ``SPOT_SOURCE_KIND``, the export pinned in code
    (``SPOT_SOURCE_EXPORTS``: name, size and sha256, also of the NAS path's name), the expected totals (spots,
    surfaces, Class A and B counts that add up to the spots, a spot diameter) and, with ``split``, the .ims of the
    split record (path, size, sha256)."""
    if not isinstance(rec, dict) or not isinstance(rec.get("changelog_entry"), int) or \
            rec["changelog_entry"] < SPOT_SOURCES_CHANGELOG_ENTRY:
        raise ValueError(f"the spot-sources record must name its changelog entry ({SPOT_SOURCES_CHANGELOG_ENTRY} or later)")
    if (rec.get("registered_split") or {}).get("roles_sha256") != roles_sha256():
        raise ValueError("the spot-sources record was made against other registered roles")
    sources = rec.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("the spot-sources record has no 'sources' list")
    tags = [s.get("tag") for s in sources]
    if len(set(tags)) != len(tags):
        raise ValueError("spot-source tags are not unique")
    by_tag = {g["tag"]: g for g in split["gonads"]} if split is not None else {}
    for s in sources:
        tag = s.get("tag")
        if tag not in PREREGISTERED_ROLES:
            raise ValueError(f"{tag}: a spot source belongs to a registered gonad")
        if PREREGISTERED_ROLES[tag] != "train" or s.get("role") != "train":
            raise ValueError(f"{tag}: a spot source is recorded for a registered train gonad only (never val or test)")
        if s.get("stem") != REGISTERED_STEMS[tag]:
            raise ValueError(f"{tag}: stem {s.get('stem')!r} is not the registered stem")
        if s.get("kind") != SPOT_SOURCE_KIND:
            raise ValueError(f"{tag}: kind {s.get('kind')!r} is not {SPOT_SOURCE_KIND!r}")
        pin = SPOT_SOURCE_EXPORTS.get(tag)
        if pin is None:
            raise ValueError(f"{tag}: no export is pinned for this gonad in the code (SPOT_SOURCE_EXPORTS)")
        ex = s.get("export") or {}
        for k in ("name", "size_bytes", "sha256"):
            if ex.get(k) != pin[k]:
                raise ValueError(f"{tag}: the export's {k} is not the pinned one")
        for k in ("path", "nas_path"):
            if ex.get(k) is not None and Path(str(ex[k]).replace("\\", "/")).name.lower() != pin["name"].lower():
                raise ValueError(f"{tag}: the export's {k} does not name the pinned file {pin['name']}")
        if not ex.get("path"):
            raise ValueError(f"{tag}: the export has no path")
        exp = s.get("expect") or {}
        if not all(isinstance(exp.get(k), int) and exp[k] > 0 for k in ("n_spots", "n_surfaces")):
            raise ValueError(f"{tag}: expect needs positive n_spots and n_surfaces")
        cc = exp.get("class_counts")
        if not isinstance(cc, dict) or not cc or sum(int(v) for v in cc.values()) != exp["n_spots"]:
            raise ValueError(f"{tag}: expect.class_counts must add up to n_spots")
        if not (isinstance(exp.get("spot_diameter_um"), list) and len(exp["spot_diameter_um"]) == 3):
            raise ValueError(f"{tag}: expect.spot_diameter_um must be [x, y, z]")
        if by_tag:
            ims = s.get("ims") or {}
            reg = by_tag[tag]["ims"]
            if (ims.get("path"), ims.get("size_bytes"), ims.get("sha256")) != (reg["path"], reg["size_bytes"], reg["sha256"]):
                raise ValueError(f"{tag}: the source's .ims is not the registered one")


def load_spot_sources(path, split: dict | None = None) -> dict:
    """Read and validate the spot-sources record (``validate_spot_sources``)."""
    rec = json.loads(Path(path).read_text(encoding="utf-8"))
    validate_spot_sources(rec, split)
    return rec


def spot_source_of(rec: dict | None, tag: str) -> dict | None:
    """The source of ``tag`` in a spot-sources record, None when it has none."""
    for s in (rec or {}).get("sources", []):
        if s["tag"] == tag:
            return s
    return None


def check_spot_source_files(rec: dict, hashes: bool = False) -> list[str]:
    """Problems with the recorded exports (missing, another size, with ``hashes`` another sha256); an empty list
    means the files on disk are the pinned ones."""
    problems = []
    for s in rec["sources"]:
        ex = s["export"]
        p = Path(ex["path"])
        if not p.is_file():
            problems.append(f"{s['tag']}: export missing: {ex['path']}")
            continue
        if p.stat().st_size != ex["size_bytes"]:
            problems.append(f"{s['tag']}: export size {p.stat().st_size} differs from the pinned {ex['size_bytes']}")
            continue
        if hashes and file_sha256(p) != ex["sha256"]:
            problems.append(f"{s['tag']}: export sha256 differs from the pinned digest")
    return problems


def spot_source_truth_usable(source: dict | None, truth_params: dict | None, truth_check: dict | None = None,
                             decision: dict | None = None) -> tuple[bool, str]:
    """Whether a registered train gonad's spot truth from its recorded earlier export (protocol changelog entry 11)
    may be scored, and why. The rules, in order, as ``addition_spot_truth_usable`` applies them to an addition:

    1. no source, or one whose export the code does not pin (``SPOT_SOURCE_EXPORTS``): not usable (the registered
       rule, region growing, stays);
    2. no truth params (no truth export from the source yet): not usable;
    3. the export: the truth params were written from this source (``spots_source`` of kind ``SPOT_SOURCE_KIND``,
       the pinned sha256), every criterion of its verification (E1 to E9) passes but the section 4.4 one (E10), and
       the frame check of entry 9 (d)(3) chose a frame. No decision overrides this;
    4. ``decision`` (of a later changelog entry): its ``usable`` and its reason;
    5. otherwise section 4.4: ``truth_check`` (this gonad's row of scripts/check_truth_assignment.py; E10's own
       agreement when no row is given) must reach ``SECTION_4_4_MIN``."""
    if source is None or SPOT_SOURCE_EXPORTS.get(source.get("tag"), {}).get("sha256") != \
            (source.get("export") or {}).get("sha256"):
        return False, "no recorded spot source (entry 11)"
    if truth_params is None:
        return False, "no truth export from the recorded spot source yet"
    src = truth_params.get("spots_source") or {}
    if src.get("source") != SPOT_SOURCE_KIND or (src.get("export") or {}).get("sha256") != source["export"]["sha256"]:
        return False, "the truth params were not written from the recorded spot source (entry 11)"
    crit = {k: v for k, v in (src.get("verification") or {}).items() if k.startswith("E")}
    if not crit:
        return False, "the recorded spot source was not verified (entry 11)"
    failed = sorted((k for k, v in crit.items() if not (v or {}).get("passes")), key=lambda k: int(k.split("_")[0][1:]))
    if any(k != EXPORT_SECTION_4_4_CRITERION for k in failed):
        return False, f"the recorded spot source failed {', '.join(k.split('_')[0] for k in failed)} of entry 11"
    fc = truth_params.get("frame_check")
    if fc is None or fc.get("chosen") is None:
        return False, "frame check of entry 9 (d)(3) missing or failed"
    agreement = None if truth_check is None else truth_check.get("agreement")
    if agreement is None:
        agreement = (crit.get(EXPORT_SECTION_4_4_CRITERION) or {}).get("agreement")
    return _decision_or_section_4_4(agreement, decision, SPOT_SOURCE_RULE)


def load_truth_check(path) -> dict[str, dict]:
    """The rows of a section 4.4 table written by scripts/check_truth_assignment.py --out, by tag, with
    ``agreement`` and the counts as numbers. ``path`` may be a list of tables (protocol changelog entry 11: the
    additions' table of entry 10 and the table of the recorded spot sources); a tag in two of them is refused."""
    import csv

    if isinstance(path, (list, tuple)):
        out: dict[str, dict] = {}
        for p in path:
            rows = load_truth_check(p)
            both = sorted(set(rows) & set(out))
            if both:
                raise ValueError(f"{p}: {both} already have a row in another section 4.4 table")
            out.update(rows)
        return out
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rec = dict(row)
            for k in ("agreement", "agreement_like_for_like"):
                try:
                    rec[k] = float(row[k]) if row.get(k) not in (None, "") else None
                except ValueError:
                    rec[k] = None
            for k in ("n_spots", "n_agree", "n_coloc_on_removed_fragment", "n_both_different"):
                try:
                    rec[k] = int(float(row[k])) if row.get(k) not in (None, "") else None
                except (KeyError, ValueError):
                    rec[k] = None
            out[str(row["tag"])] = rec
    return out


def load_spot_truth_decisions(path) -> dict:
    """The decisions file of a changelog entry on the additions' spot truth (tests/golden/rad51_additions_decisions.json,
    entry 10): ``{"changelog_entry": n, "decisions": {tag: {"usable": bool, "entry": n, "reason": str}}}``.
    Refuses a tag that is not an addition of entry 9 and an entry before 10."""
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    dec = doc.get("decisions")
    if not isinstance(dec, dict) or not dec:
        raise ValueError(f"{path}: a 'decisions' mapping of tag -> {{usable, entry, reason}} is required")
    for tag, d in dec.items():
        if tag not in ADDITION_ROLES:
            raise ValueError(f"{path}: {tag} is not an addition of protocol changelog entry 9")
        if not isinstance(d, dict) or not isinstance(d.get("usable"), bool):
            raise ValueError(f"{path}: {tag}: 'usable' must be true or false")
        if not isinstance(d.get("entry"), int) or d["entry"] < 10:
            raise ValueError(f"{path}: {tag}: 'entry' must be the changelog entry that decided it (10 or later)")
        if not str(d.get("reason") or "").strip():
            raise ValueError(f"{path}: {tag}: a reason is required")
    return doc


# ----------------------------------------------------------------------------- the manifest


def _exists_record(path, expected_size: int | None = None) -> dict:
    p = Path(path)
    exists = p.is_file()
    rec = {"path": posix(p), "exists": exists, "size_bytes": p.stat().st_size if exists else None}
    if expected_size is not None:
        rec["size_matches"] = (rec["size_bytes"] == expected_size) if exists else None
    return rec


def _addition_entry(g: dict, results_root: Path, truth_root: Path, patterns: dict,
                    truth_check: dict | None = None, decision: dict | None = None) -> dict:
    """One manifest entry for an addition of entry 9: the registered entries' fields plus ``addition``,
    ``session``, ``set``, the xlsx record, its section 4.4 row and changelog decision when given, and the reason
    for ``spot_truth_usable``."""
    stem, tag = g["stem"], g["tag"]
    rdir = results_root / tag
    outputs = {k: _exists_record(rdir / pat.format(stem=stem)) for k, pat in OUTPUT_FILES.items()}
    tdir = truth_root / tag if (truth_root / tag).is_dir() else truth_root
    truth = {k: _exists_record(tdir / pat.format(stem=stem)) for k, pat in patterns.items()}
    found = sorted(p.name for p in tdir.iterdir() if p.name.startswith(f"{stem}__imaris_")) if tdir.is_dir() else []
    params = None
    if truth.get("params", {}).get("exists"):
        try:
            params = json.loads(Path(truth["params"]["path"]).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            params = None
    usable, why = addition_spot_truth_usable(g, params, truth_check, decision)
    entry = {
        "tag": tag, "role": g["role"], "split": g["role"], "condition": g["condition"], "sex": g["sex"],
        "stem": stem, "addition": True, "session": g.get("session"), "set": g.get("set"),
        "nd2": posix(g["nd2"]["path"]), "ims": posix(g["ims"]["path"]),
        "nd2_record": _exists_record(g["nd2"]["path"], g["nd2"].get("size_bytes")),
        "ims_record": {**_exists_record(g["ims"]["path"], g["ims"].get("size_bytes")),
                       "region_growing": g["ims"].get("region_growing")},
        "results_dir": posix(rdir),
        "truth_dir": posix(truth_root / tag) if (truth_root / tag).is_dir() else posix(truth_root),
        "outputs": outputs,
        "outputs_complete": all(outputs[k]["exists"] for k in REQUIRED_OUTPUTS),
        "truth": {**truth, "found": found},
        "truth_complete": all(truth[k]["exists"] for k in REQUIRED_TRUTH if k in truth),
        "spot_truth_source": (g.get("spot_truth") or {}).get("source"),
        "spot_truth_usable": usable,
        "spot_truth_reason": why,
        "section_4_4": (None if truth_check is None else
                        {k: truth_check.get(k) for k in ("agreement", "n_spots", "n_agree",
                                                         "n_coloc_on_removed_fragment", "n_both_different")}),
        "spot_truth_decision": decision,
        "notes": list(g.get("notes", [])),
    }
    if g.get("xlsx"):
        entry["xlsx_record"] = {**_exists_record(g["xlsx"]["path"], g["xlsx"].get("size_bytes")),
                                "use": g["xlsx"].get("use")}
    return entry


def _apply_spot_source(entry: dict, source: dict, src_root: Path, patterns: dict,
                       truth_check: dict | None) -> None:
    """Point a registered entry's truth at the export of its recorded spot source (protocol changelog entry 11):
    the truth files under <src_root>/<tag>/ (the same .ims exported again, with the spots of the recorded export),
    ``spot_truth_usable`` from ``spot_source_truth_usable``, its reason, the section 4.4 row, and the source
    (``recorded_spot_source``; an addition's ``spot_truth_source`` is the name of its source, "ims" or "xlsx"); the
    registered truth folder is kept as ``registered_truth_dir`` (its nucleus truth is what the new export must
    reproduce, checked by scripts/check_truth_assignment.py --spot-sources)."""
    stem, tag = entry["stem"], entry["tag"]
    tdir = src_root / tag
    truth = {k: _exists_record(tdir / pat.format(stem=stem)) for k, pat in patterns.items()}
    found = sorted(p.name for p in tdir.iterdir() if p.name.startswith(f"{stem}__imaris_")) if tdir.is_dir() else []
    params = None
    if truth.get("params", {}).get("exists"):
        try:
            params = json.loads(Path(truth["params"]["path"]).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            params = None
    usable, why = spot_source_truth_usable(source, params, truth_check)
    entry.update({
        "registered_truth_dir": entry["truth_dir"],
        "truth_dir": posix(tdir),
        "truth": {**truth, "found": found},
        "truth_complete": all(truth[k]["exists"] for k in REQUIRED_TRUTH if k in truth),
        "recorded_spot_source": {"kind": source["kind"], "changelog_entry": SPOT_SOURCES_CHANGELOG_ENTRY,
                              "export": {k: source["export"].get(k) for k in ("name", "path", "size_bytes", "sha256")},
                              "export_record": _exists_record(source["export"]["path"], source["export"]["size_bytes"]),
                              "expect": source.get("expect")},
        "spot_truth_usable": usable,
        "spot_truth_reason": why,
        "section_4_4": (None if truth_check is None else
                        {k: truth_check.get(k) for k in ("agreement", "n_spots", "n_agree",
                                                         "n_coloc_on_removed_fragment", "n_both_different")}),
        "notes": entry["notes"] + [SPOT_SOURCE_NOTE],
    })


def build_manifest(split: dict, results_root, truth_root, name: str | None = None,
                   truth_files: dict | None = None, split_path=None, additions: dict | None = None,
                   additions_results_root=None, additions_path=None, truth_check=None,
                   spot_truth_decisions=None, spot_sources: dict | None = None, spot_sources_path=None,
                   spot_sources_truth_root=None) -> dict:
    """One entry per gonad of the split: role, condition, sex, stem, the .nd2 and .ims (exists, size
    matches the registered size), the pipeline outputs under <results_root>/<tag>/ and the truth files
    under <truth_root>, each with its path and an exists flag, plus ``outputs_complete`` (nuclei table,
    label image and spots table all present), ``truth_complete`` (spots, surfaces and masks present),
    ``spot_truth_usable`` (false when the spots were exported with region growing, or when that is
    unknown) and the split's notes. Nothing is computed from the data.

    ``additions`` (a record of ``load_additions``) appends one entry per addition after the registered ones,
    which stay exactly as they are without it: pipeline outputs under <additions_results_root>/<tag>/ (default
    the results root), truth files under <truth_root>, ``addition`` true, and ``spot_truth_usable`` from
    ``addition_spot_truth_usable`` (which reads the truth params JSON when it exists). ``truth_check`` is the
    section 4.4 table of scripts/check_truth_assignment.py (a CSV path, or a list of them) and
    ``spot_truth_decisions`` the decisions file of a changelog entry (a JSON path, ``load_spot_truth_decisions``);
    both are recorded with their sha256 under ``spot_truth_inputs``. Without them no ims-sourced addition is usable
    (entry 10 (c)).

    ``spot_sources`` (the record of ``load_spot_sources``, protocol changelog entry 11) changes only the entries of
    the registered train gonads it names: their truth files are read from <spot_sources_truth_root>/<tag>/ (the
    export from the recorded source), ``spot_truth_usable`` follows ``spot_source_truth_usable`` with this gonad's
    row of ``truth_check``, and ``recorded_spot_source`` records the source; every other registered entry stays exactly
    as without it. The record is recorded with its sha256 under ``spot_sources_file``."""
    if spot_truth_decisions is not None and additions is None:
        raise ValueError("truth_check and spot_truth_decisions apply to the additions only")
    if truth_check is not None and additions is None and spot_sources is None:
        raise ValueError("truth_check and spot_truth_decisions apply to the additions only (and truth_check to the "
                         "recorded spot sources of entry 11)")
    if spot_sources is not None and spot_sources_truth_root is None:
        raise ValueError("spot_sources needs spot_sources_truth_root, the folder of the exports from the sources")
    patterns = dict(TRUTH_FILES, **(truth_files or {}))
    results_root, truth_root = Path(results_root), Path(truth_root)
    checks = load_truth_check(truth_check) if truth_check is not None else {}
    if spot_sources is not None:
        validate_spot_sources(spot_sources, split)
    entries = []
    for g in split["gonads"]:
        stem, tag = g["stem"], g["tag"]
        rdir = results_root / tag
        outputs = {k: _exists_record(rdir / pat.format(stem=stem)) for k, pat in OUTPUT_FILES.items()}
        tdir = truth_root / tag if (truth_root / tag).is_dir() else truth_root
        truth = {k: _exists_record(tdir / pat.format(stem=stem)) for k, pat in patterns.items()}
        found = sorted(p.name for p in tdir.iterdir() if p.name.startswith(f"{stem}__imaris_")) if tdir.is_dir() else []
        rg = g["ims"].get("region_growing")
        entries.append({
            "tag": tag, "role": g["role"], "split": g["role"], "condition": g["condition"], "sex": g["sex"],
            "stem": stem,
            # plain paths under the keys germquant.benchmark.manifest.load_manifest reads; the detailed
            # records (exists, size, checksum match, region growing) sit next to them
            "nd2": posix(g["nd2"]["path"]), "ims": posix(g["ims"]["path"]),
            "nd2_record": _exists_record(g["nd2"]["path"], g["nd2"].get("size_bytes")),
            "ims_record": {**_exists_record(g["ims"]["path"], g["ims"].get("size_bytes")), "region_growing": rg},
            "results_dir": posix(rdir),
            "truth_dir": posix(truth_root / tag) if (truth_root / tag).is_dir() else posix(truth_root),
            "outputs": outputs,
            "outputs_complete": all(outputs[k]["exists"] for k in REQUIRED_OUTPUTS),
            "truth": {**truth, "found": found},
            "truth_complete": all(truth[k]["exists"] for k in REQUIRED_TRUTH if k in truth),
            "spot_truth_usable": rg is False,
            "notes": list(g.get("notes", [])),
        })
        source = spot_source_of(spot_sources, tag)
        if source is not None:
            _apply_spot_source(entries[-1], source, Path(spot_sources_truth_root), patterns, checks.get(tag))
    if additions is not None:
        validate_additions(additions)
        if (additions.get("registered_split") or {}).get("roles_sha256") != roles_sha256(
                {g["tag"]: g["role"] for g in split["gonads"]}):
            raise ValueError("the additions record was made against other roles than this split's")
        aroot = Path(additions_results_root) if additions_results_root else results_root
        dec_doc = load_spot_truth_decisions(spot_truth_decisions) if spot_truth_decisions is not None else None
        decisions = (dec_doc or {}).get("decisions", {})
        add_tags = {g["tag"] for g in additions["gonads"]}
        unknown = sorted(set(decisions) - add_tags)
        if unknown:
            raise ValueError(f"the decisions file names gonads that are not in the additions record: {unknown}")
        entries += [_addition_entry(g, aroot, truth_root, patterns, checks.get(g["tag"]), decisions.get(g["tag"]))
                    for g in additions["gonads"]]
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
    if additions is not None:
        manifest["additions_file"] = {
            "path": posix(additions_path) if additions_path else None,
            "sha256": file_sha256(additions_path) if additions_path else None,
            "changelog_entry": additions.get("changelog_entry"),
            "results_root": posix(Path(additions_results_root) if additions_results_root else results_root),
            "tags": [g["tag"] for g in additions["gonads"]]}
        manifest["spot_truth_inputs"] = {
            "truth_check": _truth_check_record(truth_check),
            "decisions": ({"path": posix(spot_truth_decisions), "sha256": file_sha256(spot_truth_decisions),
                           "changelog_entry": dec_doc.get("changelog_entry")} if dec_doc is not None else None),
            "rule": ("benchmark_split.addition_spot_truth_usable: the export's own checks, then the changelog "
                     "decision, else section 4.4 at least 0.9935 (protocol changelog entries 9 (d) and 10 (c))")}
    if spot_sources is not None:
        manifest["spot_sources_file"] = {
            "path": posix(spot_sources_path) if spot_sources_path else None,
            "sha256": file_sha256(spot_sources_path) if spot_sources_path else None,
            "changelog_entry": spot_sources.get("changelog_entry"),
            "truth_root": posix(spot_sources_truth_root),
            "tags": [s["tag"] for s in spot_sources["sources"]],
            "truth_check": _truth_check_record(truth_check),
            "rule": ("benchmark_split.spot_source_truth_usable: the export's own checks (E1 to E9 and the frame "
                     "check), then a changelog decision, else section 4.4 at least 0.9935 (protocol changelog entry "
                     "11)")}
    return manifest


def _truth_check_record(truth_check):
    """Path and sha256 of the section 4.4 table (a list of them when several were given)."""
    if truth_check is None:
        return None
    if isinstance(truth_check, (list, tuple)):
        return [{"path": posix(p), "sha256": file_sha256(p)} for p in truth_check]
    return {"path": posix(truth_check), "sha256": file_sha256(truth_check)}


def write_manifest(manifest: dict, out_root, file_name: str = "manifest.json") -> Path:
    """Write benchmarks/<name>/<file_name> (manifest.json by default) under ``out_root`` and return its path.
    A revised manifest goes next to the registered one under another file name, so the trial and final logs
    that the scripts find next to the manifest stay the same files."""
    if Path(file_name).name != file_name or not file_name.endswith(".json"):
        raise ValueError(f"{file_name!r}: a plain .json file name, no folder")
    out = Path(out_root) / manifest["name"] / file_name
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    return out


def today() -> str:
    return date.today().isoformat()
