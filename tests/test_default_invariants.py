"""The tranche 5 fence: today's default behaviour, pinned at the T5 base commit (package P0).

docs/DESIGN_system_packs.md, section 2 (rules R1 to R9) and section 11.2. Every pin is read from
tests/golden/default_invariants.json, made once at the T5 base commit on the lab workstation. A T5 package
may ADD a pin (a new recipe file arrives with one, gate G3), but an existing pin never changes, and
scripts/t5_scope_check.py fails a diff that changes or removes one. When a test here fails, today's
behaviour moved, and the fix belongs in the code, not in the JSON.

What is pinned:

* the config_hash of config/config.yaml, config_ccw77.yaml, config_n2dryice.yaml and every
  config/profiles/*.yaml, through load_config (platform independent: the text is decoded as cp1252 with
  universal newlines);
* the LF-normalised sha256 of every config/channel_maps/*.yaml and of the two Snakemake profiles under
  config/profiles/, so a CRLF checkout and an LF checkout give the same digest;
* the literal default stage hashes, stages.stage_hashes(cfg, roles, extras=FIXED, model_sha="absent"),
  for config/config.yaml with the channel roles of the two synthetic gonads, where FIXED maps every tool
  in stages._STAGE_TOOLS to "PINNED": installed versions cannot move these pins, a new payload key does;
* three default runs, each with every tool version reported as "PINNED": the synthetic and synthetic_4ch
  gonads of tests/test_golden.py (its recipe at the base is copied below, so a later edit of
  test_golden.py cannot change what the fence runs) and a single-plane (C, 1, Y, X) image read through
  the bioio branch with a fake BioImage (no bioio needed). For each run: the stage hashes it writes to
  run_manifest.json, <id>__stages.json and <id>__done.json (unlike the direct call above, these also see
  a payload key that only the pipeline's own call passes and a config value the pipeline sets before
  hashing); the key sets of run_manifest.json (and of its tool_versions), of <id>__done.json and of the
  top level of <id>__stages.json; the stage names, outcome fields and statuses of the stage record; the
  enabled stages; a few literal provenance values; the qc_flags string; the output file names; the exact
  column list of all twelve schema.TABLES (test_golden.py compares five tables and tolerates added
  columns); the sha256 of every TIF's decoded array; and the sha256 of the montage's decoded pixels;
* resume decisions (design 4.9): batch.is_done accepts the marker each default run wrote, with the
  segment hash computed as cli._batch computes it from a freshly loaded config, and also a marker with
  the base's fields;
* the legacy single-plane behaviour: today's bioio rule, spacing (1, 1, 1) with spacing_ok False even when
  the file carries its xy calibration, and a run that keeps the voxel:UNREADABLE_assumed_isotropic_1um
  flag and writes volume_um3 equal to the voxel count.

Platform rule (design 11.2): a pin that may differ between ubuntu, windows and the workstation carries
"scope": "workstation" in the JSON. It is enforced only when GERMQUANT_WORKSTATION_PINS=1, which
scripts/t5_gate.py sets on the workstation; elsewhere the test computes the digest and SKIPS with a
reason that says whether this platform agrees with the workstation pin, so a CI log (pytest -rs) records
the comparison. A pin with "scope": "all" runs everywhere. Every digest pin is its own test, so a skipped
workstation pin never hides another pin.

Adding a pin for a NEW file: run ``python tests/test_default_invariants.py OUT.json`` with PYTHONPATH=src
and copy the new entry from OUT.json into the pins file. The command writes the current values to OUT.json
and refuses to overwrite the pins file itself.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
from collections import namedtuple
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import test_end_to_end as e2e

from germquant import batch, pipeline, provenance, schema, stages
from germquant.config import load_config
from germquant.io import bioio_reader, read_image, read_image_metadata

ROOT = Path(__file__).resolve().parents[1]
PINS_FILE = ROOT / "tests" / "golden" / "default_invariants.json"
ENFORCE_WORKSTATION = os.environ.get("GERMQUANT_WORKSTATION_PINS", "") not in ("", "0")
SYNTHETIC_CHANNELS = {"synthetic": e2e.CHANNELS, "synthetic_4ch": e2e.CHANNELS4}
SINGLE_PLANE_CHANNELS = ["405", "477", "545"]
# bioio's (Z, Y, X) physical pixel sizes for a single-plane OME-TIFF: calibrated in xy, and uncalibrated
SINGLE_PLANE_PPS = ([None, 0.2, 0.2], [None, None, None])
GOLDEN_CASES = ("synthetic", "synthetic_4ch")          # the synthetic goldens of tests/test_golden.py
# The fence's default runs: case -> (image name, config tweaks). The two synthetic cases are the recipe of
# tests/test_golden.py's CASES at the T5 base (config/config.yaml, the classical segmenter, these tweaks),
# copied here so that a later edit of test_golden.py cannot change what the fence runs.
RUN_RECIPES = {
    "synthetic": ("20251105_n2_nohs_HERM_001.nd2", {}),
    "synthetic_4ch": ("20251105_n2_nohs_HERM_001.nd2", {"coloc.n_random": 20}),
    "single_plane": ("20251105_n2_nohs_HERM_001.ome.tif", {}),
}
_JUNK_SUFFIXES = (".orig", ".rej", ".bak", ".swp", "~")


def _load_pins() -> dict:
    try:
        return json.loads(PINS_FILE.read_text(encoding="utf-8"))
    except FileNotFoundError:          # only while the file is first made (``python tests/...``)
        return {}


PINS = _load_pins()


# ---------------------------------------------------------------------------------------------------
# What the pins are computed from (the generator at the bottom uses the same functions)
# ---------------------------------------------------------------------------------------------------

def lf_sha256(path: Path) -> str:
    """sha256 of the file's bytes with CRLF turned into LF: a Linux checkout, a Windows checkout with
    core.autocrlf=true (CRLF on disk) and one with core.autocrlf=false (LF on disk) agree on it."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def array_digest(a) -> dict:
    a = np.ascontiguousarray(a)
    return {"shape": list(a.shape), "dtype": str(a.dtype), "sha256": hashlib.sha256(a.tobytes()).hexdigest()}


def _rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


def _real(f: Path, root: Path) -> bool:
    """A file of the checkout under ``root``: not hidden and not in a hidden folder (.DS_Store, an editor's
    swap file, a checkpoint folder) and not a merge or editor leftover (.orig, .rej, .bak, .swp, a
    trailing ~), any of which would otherwise fail the completeness tests below."""
    return (f.is_file() and not any(part.startswith(".") for part in f.relative_to(root).parts)
            and not f.name.endswith(_JUNK_SUFFIXES))


def shipped_configs() -> list[str]:
    """config/*.yaml and config/profiles/*.yaml: every file load_config is meant to read."""
    files = sorted((ROOT / "config").glob("*.yaml")) + sorted((ROOT / "config" / "profiles").glob("*.yaml"))
    return [_rel(f) for f in files if _real(f, ROOT / "config")]


def channel_maps() -> list[str]:
    folder = ROOT / "config" / "channel_maps"
    return [_rel(f) for f in sorted(folder.glob("*.yaml")) if _real(f, folder)]


def snakemake_profiles() -> list[str]:
    """The files in the folders under config/profiles/ (the Snakemake profiles, frozen by R1 but not
    configs that load_config reads)."""
    prof = ROOT / "config" / "profiles"
    return [_rel(f) for f in sorted(prof.rglob("*")) if f.parent != prof and _real(f, prof)]


def fixed_tool_versions() -> dict[str, str]:
    return {tool: "PINNED" for tools in stages._STAGE_TOOLS.values() for tool in tools}


def pinned_tool_versions() -> dict[str, str]:
    """What provenance.tool_versions reports during a fence run: every tool it knows, as "PINNED"."""
    return {tool: "PINNED" for tool in provenance._TOOLS}


def default_stage_hashes(channel_names) -> tuple[dict, dict]:
    """(roles, stage hashes) of config/config.yaml for a gonad with these channel names, with every tool
    version replaced by "PINNED" and the model digest "absent"."""
    cfg = load_config(ROOT / "config" / "config.yaml")
    roles, _ = cfg.channel_map.resolve(list(channel_names))
    return roles, stages.stage_hashes(cfg, roles, fixed_tool_versions(), model_sha="absent")


def fence_config(case: str):
    """The run config of a fence case: config/config.yaml with the classical segmenter and the case's
    tweaks, built as tests/test_golden.py built it at the base (a later batch --resume loads the same)."""
    _, tweaks = RUN_RECIPES[case]
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    for key, value in tweaks.items():
        cfg.set(key, value)
    return cfg


def _patch_reader(mp, case: str, work: Path) -> Path:
    """Patch the reader for one fence case; returns the image path to process. synthetic and synthetic_4ch:
    the patched .nd2 reader of test_end_to_end's patched_reader / patched_reader_4ch fixtures (the same
    stack through the same two pipeline names). single_plane: a dummy .ome.tif read through the bioio
    branch, with the fake BioImage as the reader seam, so today's bioio spacing rule feeds the run."""
    image, _ = RUN_RECIPES[case]
    if case == "single_plane":
        path = work / image
        path.write_bytes(b"x")
        img = FakeBioImage(single_plane_gonad(), SINGLE_PLANE_PPS[0])
        mp.setattr(bioio_reader, "_open", lambda p: img)
        return path
    makers = {"synthetic": e2e._synthetic_stack, "synthetic_4ch": e2e._synthetic_stack_4ch}
    if case not in makers:
        raise KeyError(f"the fence has no reader for the case {case!r}")
    stack = makers[case]()
    c, z, y, x = stack.data.shape
    meta = {"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": stack.spacing, "spacing_ok": True,
            "channel_names": list(stack.channel_names), "dtype": "float32", "is_2d": False}
    mp.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    mp.setattr(pipeline, "read_stack", lambda *a, **k: stack)
    return Path(image)


def run_case(case: str, work: Path) -> dict:
    """One fence run, with every tool version reported as "PINNED", so the stage hashes it writes depend on
    their payload only. Returns the output folder, the image id, the config hash and process_image's
    result."""
    if case not in RUN_RECIPES:
        raise KeyError(f"the fence has no recipe for the case {case!r}; it runs {sorted(RUN_RECIPES)}")
    work.mkdir(parents=True, exist_ok=True)
    with pytest.MonkeyPatch.context() as mp:
        image = _patch_reader(mp, case, work)
        mp.setattr(provenance, "tool_versions", pinned_tool_versions)
        cfg = fence_config(case)
        res = pipeline.process_image(image, cfg, work / "results", prov=None)
    return {"out": Path(res["out_dir"]), "image_id": res["image_id"], "config_hash": cfg.hash, "res": res}


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def csv_columns(p: Path) -> list[str]:
    with open(p, newline="", encoding="utf-8") as fh:
        return next(csv.reader(fh))


def montage_digest(p: Path) -> dict:
    from PIL import Image

    with Image.open(p) as im:
        im.load()
        return array_digest(np.asarray(im))


def run_summary(out: Path, image_id: str) -> dict:
    """Everything the fence pins about one run's output folder."""
    import tifffile

    manifest = _json(out / "run_manifest.json")
    done = _json(out / f"{image_id}__done.json")
    record = _json(out / f"{image_id}__stages.json")
    with open(out / f"{image_id}__image_summary.csv", newline="", encoding="utf-8") as fh:
        qc_flags = next(csv.DictReader(fh))["qc_flags"]
    prefix = f"{image_id}__"
    return {
        "manifest_keys": sorted(manifest),
        "manifest_tool_version_keys": sorted(manifest["tool_versions"]),
        "done_marker_keys": sorted(done),
        "stage_record_keys": sorted(record),
        "stage_record_stages": sorted(record["stages"]),
        "stage_outcome_keys": sorted({k for oc in record["stages"].values() for k in oc}),
        "stage_status": {name: oc["status"] for name, oc in sorted(record["stages"].items())},
        "enabled_stages": list(done["enabled_stages"]),
        "provenance_values": {"model_sha256": manifest["model_sha256"], "xy_stride": manifest["xy_stride"],
                              "z_range": manifest["z_range"], "complete": record["complete"],
                              "failed_stages": done["failed_stages"]},
        "qc_flags": qc_flags,
        "files": sorted(p.name.replace(prefix, "<id>__", 1) if p.name.startswith(prefix) else p.name
                        for p in out.iterdir()),
        "table_columns": {t: csv_columns(out / f"{prefix}{t}.csv") for t in schema.TABLES},
        "tif_arrays": {p.name[len(prefix):]: array_digest(tifffile.imread(p))
                       for p in sorted(out.glob(f"{prefix}*.tif"))},
        "montage_pixels": montage_digest(out / f"{prefix}montage.png"),
        "run_stage_hashes": done["stage_hashes"],
    }


class _Dims:
    def __init__(self, sizes: dict):
        self.order = "".join(sizes)
        self.__dict__.update(sizes)


class FakeBioImage:
    """Just the surface bioio_reader uses (dims, scenes, physical_pixel_sizes, channel_names, dtype,
    reader, get_image_dask_data), so the bioio branch runs without bioio, as in tests/test_read_image.py."""

    def __init__(self, data: np.ndarray, pps):
        import dask.array as da

        self._dask = da.from_array(data, chunks=data.shape)
        self.scenes = ("Image:0",)
        self.channel_names = list(SINGLE_PLANE_CHANNELS)
        c, z, y, x = data.shape
        self.dims = _Dims({"T": 1, "C": c, "Z": z, "Y": y, "X": x})
        self.physical_pixel_sizes = namedtuple("PhysicalPixelSizes", "Z Y X")(*pps)
        self.dtype = data.dtype
        self.reader = self

    def get_image_dask_data(self, order, **kw):
        assert order == "CZYX" and kw == {"T": 0}, (order, kw)
        return self._dask


def single_plane_gonad() -> np.ndarray:
    """(C, 1, Y, X) float32 with C = DAPI, SYP, RAD-51: a row of nuclei on one plane, each with a short
    SYP stretch and one RAD-51 pixel (test_end_to_end's synthetic gonad, flattened to a single plane)."""
    y_size, x_size = 60, 200
    yy, xx = np.indices((y_size, x_size)).astype(np.float32)
    dna = np.zeros((y_size, x_size), np.float32)
    syp = np.zeros_like(dna)
    rad = np.zeros_like(dna)
    for cx in range(15, 190, 14):
        dna += np.exp(-((yy - 30) ** 2 + (xx - cx) ** 2) / (2 * 3.0 ** 2)) * 1500.0
        syp[30, cx - 3:cx + 4] = 2000.0
        rad[31, cx + 1] = 6000.0
    return np.stack([dna, syp, rad])[:, None]


def _moved(got: dict, pin: dict) -> str | None:
    """None when two stage-hash dicts agree; otherwise which stages moved (at most four named)."""
    moved = sorted(s for s in set(got) | set(pin) if got.get(s) != pin.get(s))
    if not moved:
        return None
    shown = ", ".join(f"{s} {got.get(s)} (pinned {pin.get(s)})" for s in moved[:4])
    more = f" and {len(moved) - 4} more" if len(moved) > 4 else ""
    return f"the default stage hash moved for {len(moved)} of {len(pin)} stages: {shown}{more}"


# ---------------------------------------------------------------------------------------------------
# The pins file itself
# ---------------------------------------------------------------------------------------------------

def test_pins_file_records_the_t5_base():
    assert PINS, f"{PINS_FILE} is missing: the fence has no pins"
    assert PINS["format"] == 1
    assert len(PINS["base_commit"]) == 40
    # literal anchors (P0 acceptance): config.yaml's hash and the two synthetic golden sentinels (a later
    # golden may add its own sentinel next to them)
    assert PINS["config_hashes"]["config/config.yaml"] == "84bdda5a5230"
    assert PINS["golden_sentinels"]["synthetic"] == "84bdda5a5230"
    assert PINS["golden_sentinels"]["synthetic_4ch"] == "37d49f43ac38"
    for name, sentinel in PINS["golden_sentinels"].items():
        text = (ROOT / "tests" / "golden" / name / "config_hash.txt").read_text(encoding="utf-8").strip()
        assert text == sentinel, f"tests/golden/{name}/config_hash.txt reads {text}, pinned {sentinel}"


# ---------------------------------------------------------------------------------------------------
# Frozen texts (R1): config hashes and LF-normalised digests
# ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("rel", sorted(PINS.get("config_hashes", {})))
def test_config_hash(rel):
    path = ROOT / rel
    assert path.is_file(), f"{rel} is pinned but no longer exists (rule R1: frozen files never disappear)"
    got, want = load_config(path).hash, PINS["config_hashes"][rel]
    assert got == want, (f"{rel}: config_hash {got} != pinned {want}. The text of a frozen config, its "
                         f"channel map or its base moved (rule R1); every run made with it records the old "
                         f"hash and batch --resume would reprocess it")


def test_every_shipped_config_and_profile_is_pinned():
    unpinned = sorted(set(shipped_configs()) - set(PINS["config_hashes"]))
    assert not unpinned, (f"configs or profiles without a pin: {unpinned}. A new recipe file is allowed, "
                          f"and it arrives with its own pin in {PINS_FILE.name} (gate G3)")


_TEXT_PINS = {**PINS.get("channel_maps_sha256_lf", {}), **PINS.get("snakemake_profiles_sha256_lf", {})}


@pytest.mark.parametrize("rel", sorted(_TEXT_PINS))
def test_frozen_text_digest(rel):
    path = ROOT / rel
    assert path.is_file(), f"{rel} is pinned but no longer exists (rule R1: frozen files never disappear)"
    got, want = lf_sha256(path), _TEXT_PINS[rel]
    assert got == want, (f"{rel}: LF-normalised sha256 {got[:16]}... != pinned {want[:16]}...: a frozen "
                         f"text changed, even if only a comment (rule R1)")


def test_every_channel_map_and_snakemake_profile_is_pinned():
    unpinned = sorted(set(channel_maps() + snakemake_profiles()) - set(_TEXT_PINS))
    assert not unpinned, (f"channel maps or Snakemake profiles without a pin: {unpinned}; a new file "
                          f"arrives with its own pin in {PINS_FILE.name} (gate G3)")


# ---------------------------------------------------------------------------------------------------
# Default stage hashes (risk 2: an unconditional payload key changes every default sub-hash)
# ---------------------------------------------------------------------------------------------------

def test_stage_tools():
    got = {name: list(tools) for name, tools in stages._STAGE_TOOLS.items()}
    assert got == PINS["stage_tools"], ("stages._STAGE_TOOLS changed: a new tool enters every default "
                                        "sub-hash of its stage and makes batch --resume reprocess everything")


@pytest.mark.parametrize("case", sorted(PINS.get("default_stage_hashes", {})))
def test_default_stage_hashes(case):
    pin = PINS["default_stage_hashes"][case]
    assert pin["model_sha"] == "absent"
    roles, got = default_stage_hashes(pin["channel_names"])
    assert roles == pin["roles"], f"the channel roles of the {case} gonad resolve to {roles}, pinned {pin['roles']}"
    msg = _moved(got, pin["hashes"])
    assert msg is None, (f"{case}: {msg}. The tool versions are fixed here, so the stage-hash payload itself "
                         f"changed, for example a new key on the default path (rule R2)")


# ---------------------------------------------------------------------------------------------------
# The three default runs: synthetic, synthetic_4ch and single_plane
# ---------------------------------------------------------------------------------------------------

RUN_CASES = sorted(PINS.get("runs", {}))


@pytest.fixture(scope="module")
def fence_runs(tmp_path_factory):
    runs = {}
    for case in RUN_CASES:
        run = run_case(case, tmp_path_factory.mktemp(f"fence_{case}"))
        runs[case] = {**run, "summary": run_summary(run["out"], run["image_id"])}
    return runs


@pytest.mark.parametrize("case", RUN_CASES)
def test_run_config_hash(case, fence_runs):
    run, pin = fence_runs[case], PINS["runs"][case]
    assert run["image_id"] == pin["image_id"]
    assert run["config_hash"] == pin["config_hash"]
    if case in PINS["golden_sentinels"]:
        assert pin["config_hash"] == PINS["golden_sentinels"][case]


@pytest.mark.parametrize("case", RUN_CASES)
def test_run_stage_hashes(case, fence_runs):
    """The stage hashes a default run writes, with the tool versions fixed, in all three places."""
    run, pin = fence_runs[case], PINS["runs"][case]["run_stage_hashes"]
    out, iid = run["out"], run["image_id"]
    for name in ("run_manifest.json", f"{iid}__stages.json", f"{iid}__done.json"):
        msg = _moved(_json(out / name)["stage_hashes"], pin)
        assert msg is None, (f"{case}: in {name.replace(iid, '<id>')}, {msg}. The tool versions are fixed, so "
                             f"the payload of the pipeline's own stage_hashes call changed: an argument only "
                             f"the pipeline passes, or a config value it sets before hashing (rule R2). Every "
                             f"default run would record new hashes and batch --resume could reprocess them")


@pytest.mark.parametrize("case", RUN_CASES)
def test_provenance_key_sets(case, fence_runs):
    got, pin = fence_runs[case]["summary"], PINS["runs"][case]
    for key, what in (("manifest_keys", "run_manifest.json"),
                      ("manifest_tool_version_keys", "run_manifest.json tool_versions"),
                      ("done_marker_keys", "<id>__done.json"),
                      ("stage_record_keys", "the top level of <id>__stages.json")):
        added = sorted(set(got[key]) - set(pin[key]))
        lost = sorted(set(pin[key]) - set(got[key]))
        assert not added and not lost, (f"{case}: keys of {what} changed: added {added}, removed {lost}. New "
                                        f"provenance is written only when its feature is active (rule R2)")


@pytest.mark.parametrize("case", RUN_CASES)
def test_stage_record(case, fence_runs):
    got, pin = fence_runs[case]["summary"], PINS["runs"][case]
    for key in ("stage_record_stages", "stage_outcome_keys", "stage_status", "enabled_stages",
                "provenance_values"):
        assert got[key] == pin[key], f"{case}: {key} is {got[key]}, pinned {pin[key]}"


@pytest.mark.parametrize("case", RUN_CASES)
def test_qc_flags(case, fence_runs):
    got, want = fence_runs[case]["summary"]["qc_flags"], PINS["runs"][case]["qc_flags"]
    assert got == want, (f"{case}: qc_flags {got!r} != pinned {want!r}. A qc flag on the default path moves "
                         f"the qc_flags column of every table (rule R2)")


@pytest.mark.parametrize("case", RUN_CASES)
def test_output_files(case, fence_runs):
    got, want = fence_runs[case]["summary"]["files"], PINS["runs"][case]["files"]
    added, lost = sorted(set(got) - set(want)), sorted(set(want) - set(got))
    assert not added and not lost, f"{case}: output files changed: added {added}, missing {lost}"


@pytest.mark.parametrize("case", RUN_CASES)
def test_table_columns(case, fence_runs):
    got, want = fence_runs[case]["summary"]["table_columns"], PINS["runs"][case]["table_columns"]
    assert sorted(schema.TABLES) == sorted(want), (f"schema.TABLES is now {sorted(schema.TABLES)}, pinned "
                                                   f"{sorted(want)}: a new table is written on every default run")
    bad = []
    for table, cols in want.items():
        if got[table] != cols:
            added = [c for c in got[table] if c not in cols]
            lost = [c for c in cols if c not in got[table]]
            bad.append(f"{table}: added {added}, removed {lost}" + ("" if added or lost else ", reordered"))
    assert not bad, f"{case}: table columns changed (all twelve tables are pinned): " + "; ".join(bad)


def _assert_digest(what: str, got: dict | None, pin: dict) -> None:
    """Check one digest pin. A workstation-scoped pin skips off the workstation, with a reason that says
    whether this platform agrees; each digest pin is its own test, so the skip hides nothing else."""
    assert got is not None, f"{what}: the file was not written"
    want = {k: pin[k] for k in ("shape", "dtype", "sha256")}
    same = got == want
    workstation = pin.get("scope", "all") == "workstation"
    if workstation and not ENFORCE_WORKSTATION:
        verdict = "AGREES with" if same else f"DIFFERS from (here {got['shape']} {got['sha256'][:12]})"
        pytest.skip(f"{what}: workstation-only pin, enforced when GERMQUANT_WORKSTATION_PINS=1 "
                    f"(scripts/t5_gate.py); this platform {verdict} the workstation pin {want['sha256'][:12]}")
    assert same, (f"{what}: {got['shape']} {got['dtype']} {got['sha256'][:16]}... != pinned {want['shape']} "
                  f"{want['dtype']} {want['sha256'][:16]}..."
                  + ("" if workstation else
                     " (if this differs only on another platform, design 11.2 makes the pin workstation-only)"))


TIF_PINS = [(case, name) for case in RUN_CASES for name in sorted(PINS["runs"][case].get("tif_arrays", {}))]


@pytest.mark.parametrize("case", RUN_CASES)
def test_tif_files(case, fence_runs):
    got, want = fence_runs[case]["summary"]["tif_arrays"], PINS["runs"][case]["tif_arrays"]
    assert sorted(got) == sorted(want), f"{case}: TIF files {sorted(got)}, pinned {sorted(want)}"


@pytest.mark.parametrize("case, name", TIF_PINS)
def test_tif_array(case, name, fence_runs):
    _assert_digest(f"{case} {name}", fence_runs[case]["summary"]["tif_arrays"].get(name),
                   PINS["runs"][case]["tif_arrays"][name])


@pytest.mark.parametrize("case", RUN_CASES)
def test_montage_pixels(case, fence_runs):
    _assert_digest(f"{case} montage.png decoded pixels", fence_runs[case]["summary"]["montage_pixels"],
                   PINS["runs"][case]["montage_pixels"])


# ---------------------------------------------------------------------------------------------------
# Resume decisions (design 4.9)
# ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("case", RUN_CASES)
def test_resume_accepts_the_marker_a_default_run_writes(case, fence_runs, monkeypatch):
    """batch --resume over a default run's folder, with the same config loaded afresh and the segment hash
    computed as cli._batch computes it before a batch (the tool versions fixed as in the run), must say
    "done"; otherwise every existing batch output would be reprocessed."""
    run = fence_runs[case]
    monkeypatch.setattr(provenance, "tool_versions", pinned_tool_versions)
    cfg = fence_config(case)
    model_sha = (provenance.file_sha256(pipeline.resolve_model_path(cfg))
                 if str(cfg.get("segmentation.nuclei.method", "auto")) != "classical" else "unused")
    segment_hash = stages.stage_hashes(cfg, {}, provenance.tool_versions(), model_sha,
                                       segment_backend=pipeline.segment_backend_provenance(cfg))["segment"]
    got = batch.is_done(run["out"], run["image_id"], cfg, segment_hash=segment_hash)
    assert got == (True, "done"), (f"{case}: batch.is_done says {got} for the marker this default run wrote, "
                                   f"so batch --resume would reprocess every existing default run")


def test_resume_accepts_a_marker_written_by_the_base(tmp_path):
    """A done marker with the base's fields for a default run must still mean "done" to batch --resume,
    or every existing batch output would be reprocessed."""
    run = PINS["runs"]["synthetic"]
    hashes = PINS["default_stage_hashes"]["synthetic"]["hashes"]
    marker = {"image_id": run["image_id"], "config_hash": run["config_hash"],
              "enabled_stages": run["enabled_stages"], "stage_hashes": hashes, "xy_stride": 1, "z_range": None,
              "failed_stages": [], "finished": "2026-10-11T00:00:00+00:00"}
    assert sorted(marker) == run["done_marker_keys"]
    (tmp_path / f"{run['image_id']}__done.json").write_text(json.dumps(marker), encoding="utf-8")
    cfg = load_config(ROOT / "config" / "config.yaml")
    assert batch.is_done(tmp_path, run["image_id"], cfg, segment_hash=hashes["segment"]) == (True, "done")


# ---------------------------------------------------------------------------------------------------
# The legacy single-plane behaviour
# ---------------------------------------------------------------------------------------------------

def test_single_plane_bioio_spacing_rule(tmp_path, monkeypatch):
    """Today's bioio rule for one plane: bioio reports no Z size, so the voxel size is a (1, 1, 1) guess
    with spacing_ok False, even when the file carries its xy calibration."""
    pin = PINS["single_plane"]
    p = tmp_path / "plane.ome.tif"
    p.write_bytes(b"x")
    for pps in pin["physical_pixel_sizes"]:
        img = FakeBioImage(single_plane_gonad(), pps)
        monkeypatch.setattr(bioio_reader, "_open", lambda path, img=img: img)
        st, meta = read_image(p), read_image_metadata(p)
        assert st.data.shape[1] == 1 and meta["sizes"]["Z"] == 1 and meta["is_2d"] is True
        rule = f"bioio sizes {pps} must still give spacing {pin['spacing']} with spacing_ok {pin['spacing_ok']}"
        assert list(st.spacing) == pin["spacing"] and st.spacing_ok is pin["spacing_ok"], \
            f"read_image: {rule}, got {st.spacing} with {st.spacing_ok}"
        assert list(meta["spacing"]) == pin["spacing"] and meta["spacing_ok"] is pin["spacing_ok"], \
            f"read_image_metadata: {rule}, got {meta['spacing']} with {meta['spacing_ok']}"


def test_single_plane_run_keeps_the_legacy_flag_and_volumes(fence_runs):
    """The single-plane fence run (a calibrated OME-TIFF through process_image's bioio branch) flags the
    voxel size and writes every volume as a voxel count in um3. Its exact qc_flags string, key sets, file
    names, columns, stage hashes and arrays are pinned with the other runs above."""
    pin = PINS["single_plane"]
    run = fence_runs["single_plane"]
    flags = run["res"]["qc_flags"]
    assert pin["flag"] in flags, (f"a single-plane image lost the legacy {pin['flag']} flag (flags: {flags}); "
                                  f"planar behaviour is opt-in (design 4.1)")
    summary = pd.read_csv(run["out"] / f"{run['image_id']}__image_summary.csv")
    assert pin["flag"] in str(summary["qc_flags"].iloc[0]).split(";")
    nuclei = pd.read_csv(run["out"] / f"{run['image_id']}__nuclei.csv")
    assert len(nuclei) > 0
    assert (nuclei["volume_um3"] == nuclei["n_voxels"]).all(), ("a single-plane volume_um3 is no longer the voxel "
                                                                "count of today's (1, 1, 1) guess")
    for col in ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um"):
        assert (nuclei[col] == 1.0).all(), f"{col} is no longer today's 1.0 um guess for a single plane"


# ---------------------------------------------------------------------------------------------------
# Generator: writes the current values to a file of your choice (never to the pins file)
# ---------------------------------------------------------------------------------------------------

ABOUT = ("The tranche 5 fence (docs/DESIGN_system_packs.md 11.2), read by tests/test_default_invariants.py. "
         "Made once at base_commit on the lab workstation. Pins may be added for new files; existing pins "
         "never change (scripts/t5_scope_check.py fails a diff that changes or removes one). scope "
         "'workstation' marks a pin that is enforced only on the workstation (GERMQUANT_WORKSTATION_PINS=1, "
         "set by scripts/t5_gate.py).")


def _environment() -> dict:
    env = {"platform": platform.platform(), "python": platform.python_version()}
    for dist in ("numpy", "scipy", "scikit-image", "pandas", "tifffile", "matplotlib", "pillow", "spotmax"):
        try:
            env[dist] = metadata.version(dist)
        except metadata.PackageNotFoundError:
            env[dist] = "not-installed"
    return env


def current_pins(work: Path, base_commit: str) -> dict:
    pins: dict = {"format": 1, "about": ABOUT, "base_commit": base_commit, "made_on": _environment()}
    pins["config_hashes"] = {rel: load_config(ROOT / rel).hash for rel in shipped_configs()}
    pins["golden_sentinels"] = {name: (ROOT / "tests" / "golden" / name / "config_hash.txt")
                                .read_text(encoding="utf-8").strip() for name in GOLDEN_CASES}
    pins["channel_maps_sha256_lf"] = {rel: lf_sha256(ROOT / rel) for rel in channel_maps()}
    pins["snakemake_profiles_sha256_lf"] = {rel: lf_sha256(ROOT / rel) for rel in snakemake_profiles()}
    pins["stage_tools"] = {name: list(tools) for name, tools in stages._STAGE_TOOLS.items()}
    pins["default_stage_hashes"] = {}
    for name, names in SYNTHETIC_CHANNELS.items():
        roles, hashes = default_stage_hashes(names)
        pins["default_stage_hashes"][name] = {"channel_names": list(names), "roles": roles,
                                              "model_sha": "absent", "hashes": hashes}
    pins["runs"] = {}
    for case in RUN_RECIPES:
        run = run_case(case, work / case)
        summary = run_summary(run["out"], run["image_id"])
        summary["tif_arrays"] = {k: {**v, "scope": "all"} for k, v in summary["tif_arrays"].items()}
        summary["montage_pixels"] = {**summary["montage_pixels"], "scope": "workstation"}
        pins["runs"][case] = {"image_id": run["image_id"], "config_hash": run["config_hash"], **summary}
    pins["single_plane"] = {"flag": "voxel:UNREADABLE_assumed_isotropic_1um",
                            "physical_pixel_sizes": [list(p) for p in SINGLE_PLANE_PPS],
                            "spacing": [1.0, 1.0, 1.0], "spacing_ok": False}
    return pins


def dump_pins(obj, indent: int = 0) -> str:
    """JSON with one key per line and short lists of scalars kept on one line."""
    pad, inner = " " * indent, " " * (indent + 1)
    if isinstance(obj, dict) and obj:
        items = [f"{inner}{json.dumps(k)}: {dump_pins(v, indent + 1)}" for k, v in obj.items()]
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    return json.dumps(obj)


def main(argv: list[str]) -> int:
    if len(argv) != 1:
        print("usage: python tests/test_default_invariants.py OUT.json  (writes the current values to OUT.json)")
        return 2
    out = Path(argv[0]).resolve()
    if out == PINS_FILE.resolve():
        print(f"refusing to write {PINS_FILE}: existing pins never change; copy a NEW entry by hand")
        return 2
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True,
                          check=True).stdout.strip()
    with tempfile.TemporaryDirectory() as work:
        pins = current_pins(Path(work), PINS.get("base_commit", head))
    out.write_text(dump_pins(pins) + "\n", encoding="utf-8")
    print(f"current values written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
