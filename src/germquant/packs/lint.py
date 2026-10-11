"""`germquant pack lint DIR [--previous ZIP]` (docs/DESIGN_system_packs.md, section 3.7).

Every pack: ASCII text files only (so no em or en dashes); no config/config.yaml inside the pack (it would
shadow a bare ``extends: config.yaml``); each recipe's expected_config_hash equals its computed hash; every
tunable is a key the recipe sets, with a scalar value (what packs.validate_pack refuses); SHA256SUMS is
current; pack: lines are versioned; no absolute or machine path in a recipe, a channel map, a card or
pack.yaml; every card agrees with its model spec (models_registry.validate_card); a root: checkout pack's
engine-file digests are current; a local example image exists and stays under 2 MB (a real image is fetched by
`germquant pack fetch-example`, never shipped), and a real example's base_url is https.

Outside recipes (root: pack) also:

* set every worm stage switch explicitly (germline, envelope, audit, axis, staging, sc.trace, granule, coloc,
  partition, granule_tail);
* set the neutral metadata.filename_regex ``(?!x)x``, which matches nothing, unless pack.yaml declares the
  pack's own regex under ``filename_regex``;
* set segmentation.nuclei.method explicitly to cellpose (with any backend) or classical, never auto, so a
  missing Cellpose fails loudly instead of falling back to the watershed;
* reach no config/config.yaml through their extends chain (they extend the pack's frozen pack_base.yaml);
* switch the batch framing check off through qc.germline_outlier_factor (at least 1e6), set qc.min_nuclei,
  and set expected_from_germ_cell: false on every spots instance;
* set only keys the engine reads: a key config/config.yaml sets (read from the engine checkout) or a key on
  ALLOW_KEYS / ALLOW_PREFIXES (code-default keys config.yaml does not carry, and the keys tranche 5 adds).
  DENY_KEYS wins over config.yaml: keys it carries that nothing reads (spots.channel, the spots stage always
  counts the foci role, pipeline.py). A misspelt key fails.

The C. elegans pack (root: checkout) is exempt from the outside rules: its recipes are the frozen engine
files themselves. ``--previous`` checks the version bump against the previous release (semver rules below).
"""
from __future__ import annotations

import json
import re
import tempfile
import zipfile
from pathlib import Path

import yaml

from .format import (
    MANIFEST,
    TEXT_SUFFIXES,
    Pack,
    PackFormatError,
    load_pack_folder,
    pack_files,
    sums_sha256,
    verify_sums,
)

WORM_SWITCHES = ("germline.enabled", "envelope.enabled", "audit.enabled", "axis.enabled", "staging.enabled",
                 "sc.trace.enabled", "granule.enabled", "coloc.enabled", "partition.enabled", "granule_tail.enabled")
NEUTRAL_REGEX = "(?!x)x"
FRAMING_OFF_MIN = 1.0e6
EXAMPLE_MAX_BYTES = 2_000_000          # a pack's own (synthetic) example; real images are fetched (design 3.2, 6.9)
EN_DASH, EM_DASH = chr(0x2013), chr(0x2014)      # refused in pack text (kept out of this source)
# config.yaml carries these, but nothing reads them (the spots stage always counts the foci role)
DENY_KEYS = frozenset({"spots.channel", "project", "axis.method", "output.dir"})
# keys the engine reads with a code default that config/config.yaml does not carry
ALLOW_KEYS = frozenset({
    "acquisition.read_exposures", "signal_qc.enabled", "germline.enabled", "envelope.enabled", "audit.enabled",
    "axis.enabled", "staging.enabled", "sc.trace.enabled", "granule.enabled", "partition.enabled",
    "granule_tail.enabled",
    "segmentation.nuclei.backend", "segmentation.nuclei.device", "segmentation.nuclei.channel",
    "spots.backend", "spots.instances",
    "envelope.band_um", "envelope.smooth_um", "envelope.volume_ratio_min", "envelope.volume_ratio_max",
    "envelope.ring_smooth_um", "envelope.ring_shell_um", "envelope.ring_ratio_max", "envelope.shell_over_thr_max",
    "envelope.territory_dilate_um", "envelope.cyto_um", "audit.missed_cov_max",
    "granule.method", "granule.region",
    "staging.traces_file", "staging.centroid", "staging.off_axis_um", "staging.adaptive_cutoff",
    "sc.intensity_percentile", "sc.min_fragment_length_um",
    "partition.bg_percentile", "partition.zshift_planes", "partition.whole_min_cyto_voxels",
    "partition.whole_min_granule_objects", "partition.min_syp_snr",
    "signal_qc.stride", "signal_qc.tissue_dilation_um", "signal_qc.ring_um", "signal_qc.min_bg_voxels",
    "signal_qc.bit_depth", "signal_qc.flags.enabled", "signal_qc.min_snr", "signal_qc.max_sat_frac",
    "signal_qc.crosstalk.max_r",
    "output.write_envelope_labels", "qc.exclusions_file",
})
# keys tranche 5 adds (docs/DESIGN_system_packs.md, sections 3.10, 3.12 and 4.2), read by the packages that add them
T5_KEYS = frozenset({
    "io.planar", "io.voxel_size_um", "io.input_suffixes",
    "segmentation.nuclei.project", "segmentation.nuclei.project_plane", "segmentation.nuclei.min_area_um2",
    "segmentation.nuclei.cellpose.object_diameter_um", "segmentation.nuclei.labels_from",
    "segmentation.nuclei.require_channel",
    "spots.spot_radius_z_um", "spots.remove_hot_pixels", "spots.sharpen", "qc.min_nuclei",
})
# sections whose inner keys belong to a backend or are data (a mapping read as a whole)
ALLOW_PREFIXES = ("metadata.defaults", "acquisition.expected", "spots.log", "spots.net3d",
                  "segmentation.nuclei.nuclei3d", "granule.tophat", "timelapse")


def engine_keys(config_yaml: Path) -> set[str]:
    """The dotted leaf keys config/config.yaml sets (lists are leaves)."""
    data = yaml.safe_load(Path(config_yaml).read_text(encoding="cp1252")) or {}
    return set(_flatten(data))


def _flatten(d, prefix: str = "") -> list[str]:
    out = []
    for k, v in (d or {}).items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict) and v and not any(key == p or key.startswith(p + ".") for p in ALLOW_PREFIXES):
            out.extend(_flatten(v, key))
        else:
            out.append(key)
    return out


def key_readable(key: str, known: set[str]) -> bool:
    if key in DENY_KEYS:
        return False
    if key in known or key in ALLOW_KEYS or key in T5_KEYS:
        return True
    return any(key == p or key.startswith(p + ".") for p in ALLOW_PREFIXES)


def _get(d: dict, dotted: str, default=None):
    node = d
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node


def _string_leaves(obj, path: str = ""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _string_leaves(v, f"{path}.{k}" if path else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _string_leaves(v, f"{path}[{i}]")
    elif isinstance(obj, str):
        yield path, obj


def _machine_path(value: str) -> bool:
    from ..models_registry import _is_machine_path

    v = value.strip()
    if "://" in v:
        return False
    return _is_machine_path(v) or v.startswith("~/") or v.startswith("~\\")


def _chain_files(recipe_path: Path) -> list[Path]:
    """The files of a recipe's extends chain, the recipe first (resolved as load_config resolves them)."""
    from ..config import _base_dir_for, _resolve_extends_ex
    from ..fsutil import read_config_text

    files, seen = [], set()
    cur = Path(recipe_path)
    while cur is not None and cur.resolve() not in seen:
        seen.add(cur.resolve())
        files.append(cur)
        data = yaml.safe_load(read_config_text(cur)) or {}
        ext = data.get("extends")
        if not ext:
            break
        cur, _ = _resolve_extends_ex(cur, _base_dir_for(cur), str(ext))
    return files


def _inside(path: Path, folder: Path) -> bool:
    try:
        path.resolve().relative_to(folder.resolve())
        return True
    except ValueError:
        return False


def lint_pack(folder: str | Path, previous: str | Path | None = None) -> dict:
    """Run every lint rule on the pack in ``folder``; ``previous`` is the previous release's zip. Returns
    ``{"pack", "problems", "warnings"}``; the pack is clean when ``problems`` is empty."""
    from ..config import builtin_roots, load_config
    from ..models_registry import validate_card
    from .discovery import engine_root_for, recipe_profile_path
    from .status import engine_files_problems
    from .validate import tunable_problems

    folder = Path(folder)
    problems: list[str] = []
    warnings: list[str] = []
    try:
        pack = load_pack_folder(folder)
    except PackFormatError as e:
        return {"pack": str(folder), "problems": [str(e)], "warnings": []}
    problems += [f"manifest: {p}" for p in pack.problems]

    for rel in pack_files(folder):                                  # ASCII, no dashes
        if Path(rel).suffix.lower() not in TEXT_SUFFIXES and rel != "SHA256SUMS":
            continue
        raw = (folder / rel).read_bytes()
        try:
            raw.decode("ascii")
        except UnicodeDecodeError as e:
            text = raw.decode("utf-8", "replace")
            dash = " (an en or em dash)" if (EN_DASH in text or EM_DASH in text) else ""
            problems.append(f"{rel}: byte {e.start} is not ASCII{dash}; pack files are ASCII only")
    if (folder / "config" / "config.yaml").exists():
        problems.append("config/config.yaml inside the pack would shadow a bare `extends: config.yaml`")
    problems += [f"SHA256SUMS: {p}" for p in verify_sums(folder)]
    problems += engine_files_problems(pack)

    for path_rel, value in _string_leaves(pack.data):              # pack.yaml
        if _machine_path(value):
            problems.append(f"pack.yaml {path_rel}: a machine path {value!r}")
    for ex in pack.data.get("examples") or []:                       # a local example stays tiny (design 3.2)
        if isinstance(ex, dict) and ex.get("sources") and not str(ex.get("base_url", "")).startswith("https://"):
            problems.append(f"example {ex.get('id')!r}: base_url {ex.get('base_url')!r} is not https (plain http to "
                            f"this machine is for tests only)")
        if not isinstance(ex, dict) or not ex.get("image"):
            continue
        img = folder / str(ex["image"])
        if not img.is_file():
            problems.append(f"example {ex.get('id')!r}: {ex['image']} is missing")
        elif img.stat().st_size >= EXAMPLE_MAX_BYTES:
            problems.append(f"example {ex.get('id')!r}: {ex['image']} is {img.stat().st_size:,} bytes; a pack's own "
                            f"example stays under {EXAMPLE_MAX_BYTES:,} (real images are fetched, never shipped)")
    for spec in pack.models:                                         # cards
        if not spec.card:
            warnings.append(f"model {spec.name}: no card")
            continue
        cp = folder / spec.card
        if not cp.is_file():
            problems.append(f"model {spec.name}: card {spec.card} is missing")
            continue
        try:
            card = json.loads(cp.read_text(encoding="utf-8"))
        except ValueError as e:
            problems.append(f"{spec.card}: not valid JSON ({e})")
            continue
        res = validate_card(card, spec.as_dict(), pack_written=True)
        problems += [f"{spec.card}: {p}" for p in res["problems"]]
        problems += [f"{spec.card}: {k} holds a machine path" for k in res["non_portable"]]
        for path_rel, value in _string_leaves(card):
            if _machine_path(value) and path_rel not in ("model_file", "base_model.path"):
                problems.append(f"{spec.card} {path_rel}: a machine path {value!r}")

    own_yaml = [folder / rel for rel in pack_files(folder) if rel.startswith("config/")
                and Path(rel).suffix.lower() in (".yaml", ".yml")]
    for f in own_yaml:                                               # recipes and maps of the pack itself
        try:
            data = yaml.safe_load(f.read_text(encoding="ascii", errors="replace")) or {}
        except yaml.YAMLError as e:
            problems.append(f"{f.relative_to(folder).as_posix()}: not valid YAML ({e})")
            continue
        rel = f.relative_to(folder).as_posix()
        ext = data.get("extends") if isinstance(data, dict) else None
        if isinstance(ext, str) and ext.startswith("pack:") and "@" not in ext.split("/", 1)[0]:
            problems.append(f"{rel}: `extends: {ext}` must name a pack version (pack:<id>@<version>/<recipe>)")
        for path_rel, value in _string_leaves(data):
            if _machine_path(value):
                problems.append(f"{rel} {path_rel}: a machine path {value!r}")

    known: set[str] | None = None
    if pack.root == "pack":
        roots = builtin_roots()
        cfg_yaml = roots[0] / "config" / "config.yaml" if roots else None
        if cfg_yaml is not None and cfg_yaml.is_file():
            known = engine_keys(cfg_yaml)
        else:
            warnings.append("no engine checkout found: the key rule could not read config/config.yaml")
    for r in pack.recipes:
        try:
            path = recipe_profile_path(pack, r)
            cfg = load_config(path)
        except Exception as e:  # noqa: BLE001 - reported as a problem
            problems.append(f"recipe {r.id}: does not load ({type(e).__name__}: {e})")
            continue
        if cfg.hash != r.expected_config_hash:
            problems.append(f"recipe {r.id}: expected_config_hash {r.expected_config_hash} but it hashes to {cfg.hash}")
        merged = cfg.as_dict()
        problems += tunable_problems(r, merged)                     # what validate_pack refuses
        if pack.root == "checkout":
            continue
        problems += [f"recipe {r.id}: {p}" for p in _outside_rules(pack, r, path, merged, known)]
    if previous is not None:
        problems += _previous_rules(pack, Path(previous))
    try:
        if pack.root == "checkout":
            engine_root_for(pack)
    except FileNotFoundError as e:
        warnings.append(str(e))
    return {"pack": pack.ref, "problems": problems, "warnings": warnings}


def _outside_rules(pack: Pack, recipe, path: Path, merged: dict, known: set[str] | None) -> list[str]:
    out: list[str] = []
    try:
        chain = _chain_files(path)
    except Exception as e:  # noqa: BLE001
        return [f"extends chain does not resolve ({type(e).__name__}: {e})"]
    for f in chain:
        if f.name == "config.yaml" and f.parent.name == "config":
            out.append(f"the extends chain reaches {f}: outside recipes extend the pack's frozen pack_base.yaml, "
                       f"never config/config.yaml")
    own: dict = {}
    for f in chain:
        if not _inside(f, pack.folder):
            continue
        data = yaml.safe_load(f.read_text(encoding="ascii", errors="replace")) or {}
        data.pop("extends", None)
        for key in _flatten(data):
            own.setdefault(key, f)
            if known is not None and not key_readable(key, known):
                why = "nothing reads it" if key in DENY_KEYS else "the engine does not read it (a typo?)"
                out.append(f"{f.relative_to(pack.folder).as_posix()}: key {key}: {why}")
    for sw in WORM_SWITCHES:
        if sw not in own:
            out.append(f"worm stage switch {sw} is not set explicitly")
    want_regex = str(pack.data.get("filename_regex", NEUTRAL_REGEX))
    if _get(merged, "metadata.filename_regex") != want_regex:
        out.append(f"metadata.filename_regex must be {want_regex!r} (the neutral regex unless pack.yaml declares "
                   f"filename_regex); the worm regex misreads other file names")
    method = _get(merged, "segmentation.nuclei.method")
    if "segmentation.nuclei.method" not in own or method not in ("cellpose", "classical"):
        out.append(f"segmentation.nuclei.method must be set to cellpose (any backend) or classical, never auto "
                   f"(got {method!r})")
    factor = _get(merged, "qc.germline_outlier_factor")
    if not isinstance(factor, (int, float)) or isinstance(factor, bool) or factor < FRAMING_OFF_MIN:
        out.append(f"qc.germline_outlier_factor must switch the batch framing check off (>= {FRAMING_OFF_MIN:g}, for "
                   f"example 1.0e+9), got {factor!r}")
    mn = _get(merged, "qc.min_nuclei")
    if isinstance(mn, bool) or not isinstance(mn, int) or mn < 0:
        out.append(f"qc.min_nuclei must be set to a count (0 or more), got {mn!r}")
    for i, inst in enumerate(_get(merged, "spots.instances", []) or []):
        if not isinstance(inst, dict) or inst.get("expected_from_germ_cell", True) is not False:
            out.append(f"spots.instances[{i}] must set expected_from_germ_cell: false (a worm-only column otherwise)")
    return out


def _recipe_hashes(pack: Pack) -> dict[str, str]:
    from ..config import load_config
    from .discovery import recipe_profile_path

    out = {}
    for r in pack.recipes:
        try:
            out[r.id] = load_config(recipe_profile_path(pack, r)).hash
        except Exception:  # noqa: BLE001 - an unloadable recipe counts as changed
            out[r.id] = "unloadable"
    return out


def _semver(v: str) -> tuple[int, int, int]:
    m = re.match(r"^(\d+)\.(\d+)\.(\d+)", v)
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else (0, 0, 0)


def required_bump(old: Pack, new: Pack) -> tuple[str | None, list[str]]:
    """(``major`` | ``minor`` | ``patch`` | None, reasons): the smallest bump the changes need. Before 1.0.0
    any change to hashed text or a model digest needs a minor bump; from 1.0.0 it needs a major one. A new
    recipe or a widened engine range needs a minor bump; anything else (card, docs, metadata) a patch."""
    reasons: list[str] = []
    level = None
    oh, nh = _recipe_hashes(old), _recipe_hashes(new)
    om = {m.name: m.sha256 for m in old.models}
    nm = {m.name: m.sha256 for m in new.models}
    hashed = [rid for rid in oh if rid in nh and oh[rid] != nh[rid]]
    removed = [rid for rid in oh if rid not in nh]
    digests = [n for n in om if n in nm and om[n] != nm[n]]
    if hashed or removed or digests:
        reasons += [f"recipe {rid} hashes differently" for rid in hashed]
        reasons += [f"recipe {rid} was removed" for rid in removed]
        reasons += [f"model {n} has another digest" for n in digests]
        level = "major" if _semver(old.version)[0] >= 1 else "minor"
    added = [rid for rid in nh if rid not in oh]
    widened = (old.data.get("requires") or {}).get("germquant") != (new.data.get("requires") or {}).get("germquant")
    if added or widened:
        reasons += [f"new recipe {rid}" for rid in added] + (["the engine range changed"] if widened else [])
        level = level or "minor"
    if level is None and sums_sha256(old.folder) != sums_sha256(new.folder):
        reasons.append("card, docs or metadata changed")
        level = "patch"
    return level, reasons


def _previous_rules(pack: Pack, previous: Path) -> list[str]:
    from .install import PackInstallError, _pack_folder_in, extract_zip

    with tempfile.TemporaryDirectory() as tmp:
        try:
            if previous.is_dir():
                folder = _pack_folder_in(previous)
            elif zipfile.is_zipfile(previous):
                extract_zip(previous, Path(tmp))
                folder = _pack_folder_in(Path(tmp))
            else:
                return [f"--previous {previous}: neither a pack folder nor a zip"]
            old = load_pack_folder(folder)
        except (PackInstallError, PackFormatError, OSError) as e:
            return [f"--previous {previous}: {e}"]
        if old.id != pack.id:
            return [f"--previous is pack {old.id}, not {pack.id}"]
        level, reasons = required_bump(old, pack)
        o, n = _semver(old.version), _semver(pack.version)
    if level is None:
        return [] if n >= o else [f"version {pack.version} is lower than the previous {old.version}"]
    if n <= o:
        return [f"version {pack.version} must be above the previous {old.version}: " + "; ".join(reasons)]
    ok = {"major": n[0] > o[0],
          "minor": n[0] > o[0] or (n[0] == o[0] and n[1] > o[1]),
          "patch": n > o}[level]
    if not ok:
        return [f"the changes need a {level} version bump over {old.version} (got {pack.version}): " + "; ".join(reasons)]
    return []


def manifest_text(folder: Path) -> str:
    return (Path(folder) / MANIFEST).read_text(encoding="ascii", errors="replace")
