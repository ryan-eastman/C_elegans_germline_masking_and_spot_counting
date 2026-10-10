"""The generic profile (config/profiles/generic.yaml with config/channel_maps/generic_nuclei_puncta.yaml):
samples that are not a worm germline. Covers the map's roles, the profile's hash rule, its stage
switches, a two-channel (nuclear stain + puncta) synthetic stack end to end on the CPU, and the
in_germline contract when the germline stage is off. docs/GENERIC_PROFILE.md."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel worm gonad)

from germquant import pipeline, schema
from germquant.config import load_config
from germquant.io.channel_map import ChannelMap
from germquant.io.nd2_reader import Stack
from germquant.stages import stage_enabled

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "config" / "profiles" / "generic.yaml"
CHANNEL_MAP = ROOT / "config" / "channel_maps" / "generic_nuclei_puncta.yaml"
SPACING = (0.4, 0.2, 0.2)
CHANNELS2 = ["DAPI", "488"]          # a nuclear stain and a puncta channel named by laser line
WORM_ONLY = ("germline", "axis", "staging", "envelope", "audit", "sc_trace", "partition", "granule_tail")
ALWAYS = ("read", "segment", "measure", "spots", "qc", "render", "write")


def _synthetic_stack_2ch():
    """Nuclei along a line, one planted punctum each: built like test_end_to_end's gonad but with no
    SYP channel and no off-gonad junk (there is no germline to isolate in a dish)."""
    Z, Y, X = 20, 80, 320
    dapi = np.zeros((Z, Y, X), np.float32)
    foci = np.zeros((Z, Y, X), np.float32)
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    for cx in range(20, 300, 17):                              # 17 nuclei along x at y=40
        r2 = (((zz - 10) * SPACING[0]) ** 2 + ((yy - 40) * SPACING[1]) ** 2
              + ((xx - cx) * SPACING[2]) ** 2)
        dapi += np.exp(-r2 / (2 * 1.6 ** 2)) * 1500.0
        foci[10, 42, min(cx + 1, X - 1)] = 6000.0              # one punctum per nucleus
    return Stack(data=np.stack([dapi, foci], axis=0), spacing=SPACING, channel_names=CHANNELS2,
                 path="x.nd2", spacing_ok=True)


@pytest.fixture
def patched_reader_2ch(monkeypatch):
    stack = _synthetic_stack_2ch()
    meta = {"sizes": {"C": 2, "Z": 20, "Y": 80, "X": 320}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS2, "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)


def test_generic_channel_map_roles():
    cm = ChannelMap.from_yaml(CHANNEL_MAP)
    assert list(cm.roles) == ["dna", "foci", "granule"]
    assert cm.roles["dna"].required
    assert not cm.roles["foci"].required and not cm.roles["granule"].required
    assert not cm.has_role("central_element")
    # descriptive names resolve by name, without a flag
    roles, flags = cm.resolve(["DNA", "RAD51 foci", "stress granules"])
    assert roles == {"dna": 0, "foci": 1, "granule": 2} and flags == []
    # laser-line names: the position fallback, flagged per role that fell back
    roles, flags = cm.resolve(CHANNELS2)
    assert roles == {"dna": 0, "foci": 1, "granule": None}
    assert flags == ["channel:foci:name_unmatched_used_index_1"]
    # a nuclear stain alone: the optional roles stay unresolved without any flag
    roles, flags = cm.resolve(["Hoechst"])
    assert roles == {"dna": 0, "foci": None, "granule": None} and flags == []
    # one unnamed channel still becomes dna by position (flagged); no channel at all is a hard miss
    roles, flags = cm.resolve(["GFP"])
    assert roles["dna"] == 0 and "channel:dna:name_unmatched_used_index_0" in flags
    assert "channel:dna:MISSING_REQUIRED" in cm.resolve([])[1]
    assert cm.marker("foci") == "foci" and cm.marker("dna") == "DNA"


def test_generic_profile_loads_and_keeps_the_base_hash():
    sentinel = (ROOT / "tests" / "golden" / "synthetic" / "config_hash.txt").read_text().strip()
    base = load_config(ROOT / "config" / "config.yaml")
    assert base.hash == sentinel
    gen = load_config(PROFILE)
    assert gen.hash != base.hash                                  # its own recipe hash ...
    assert gen.hash == load_config(PROFILE).hash                  # ... stable across loads
    assert load_config(ROOT / "config" / "config.yaml").hash == sentinel   # the base is untouched
    assert "extends" not in gen.as_dict() and gen.base_dir == ROOT
    assert gen.get("io.channel_map") == "config/channel_maps/generic_nuclei_puncta.yaml"
    assert set(gen.channel_map.roles) == {"dna", "foci", "granule"}
    assert gen.get("segmentation.nuclei.cellpose_model") == "cpsam"
    assert gen.get("io.exclude_patterns") == []
    for name in WORM_ONLY:
        assert stage_enabled(gen, name) is False, name
    for name in ALWAYS:
        assert stage_enabled(gen, name) is True, name
    # granule / coloc are left to their channel-role gate, not switched off by the profile
    assert stage_enabled(gen, "coloc") is True and stage_enabled(gen, "granule") is True
    # the spot recipe is inherited from the base (one recipe, not a fork of it)
    assert gen.get("spots.effect_size_min") == base.get("spots.effect_size_min")
    assert gen.get("spots.thresholding_method") == base.get("spots.thresholding_method")
    assert gen.get("spots.spot_radius_um") == base.get("spots.spot_radius_um")


def test_study_profile_chains_on_top_of_generic(tmp_path):
    """docs/GENERIC_PROFILE.md: a study overlay extends generic.yaml (which extends config.yaml). Out of
    the repository it names generic.yaml by absolute path and inherits the generic channel map as an
    absolute path; its hash differs from the generic profile's."""
    gen = load_config(PROFILE)
    study = tmp_path / "hela.yaml"
    study.write_text(f"extends: {PROFILE.as_posix()}\nspots:\n  spot_radius_um: 0.25\n", encoding="utf-8")
    cfg = load_config(study)
    assert cfg.base_dir == tmp_path
    cm = Path(str(cfg.get("io.channel_map")))
    assert cm.is_absolute() and cm.is_file() and cm.name == CHANNEL_MAP.name
    assert set(cfg.channel_map.roles) == {"dna", "foci", "granule"}
    assert cfg.get("spots.spot_radius_um") == 0.25
    assert cfg.get("spots.effect_size_min") == gen.get("spots.effect_size_min")
    assert cfg.get("segmentation.nuclei.cellpose_model") == "cpsam"
    for name in WORM_ONLY:
        assert stage_enabled(cfg, name) is False, name
    assert cfg.hash != gen.hash and cfg.hash == load_config(study).hash
    # inside config/profiles/ a bare name is found next to the study file
    d = tmp_path / "config" / "profiles"
    d.mkdir(parents=True)
    (d / "generic.yaml").write_text(PROFILE.read_text(encoding="utf-8").replace(
        "extends: config.yaml", f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}"), encoding="utf-8")
    (tmp_path / "config" / "channel_maps").mkdir()
    (tmp_path / "config" / "channel_maps" / CHANNEL_MAP.name).write_text(CHANNEL_MAP.read_text(encoding="utf-8"),
                                                                        encoding="utf-8")
    (d / "study.yaml").write_text("extends: generic.yaml\nspots:\n  spot_radius_um: 0.2\n", encoding="utf-8")
    cfg2 = load_config(d / "study.yaml")
    assert cfg2.get("spots.spot_radius_um") == 0.2 and stage_enabled(cfg2, "germline") is False
    assert set(cfg2.channel_map.roles) == {"dna", "foci", "granule"}


def test_generic_profile_runs_a_two_channel_stack_end_to_end(tmp_path, patched_reader_2ch):
    cfg = load_config(PROFILE)
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"     # CPU, deterministic (as every e2e test)
    out = tmp_path / "generic"
    res = pipeline.process_image("plate1_wellA01_site3.nd2", cfg, out, prov=None)

    assert res["n_nuclei"] >= 10
    assert not any(":FAILED_" in f for f in res["qc_flags"]), res["qc_flags"]
    assert not [s for s, o in res["stages"].items() if o["status"] == "failed"], res["stages"]
    st = res["stages"]
    for name in ALWAYS:
        assert st[name]["status"] == "ran", name
    for name in WORM_ONLY:
        assert (st[name]["status"], st[name]["reason"]) == ("skipped", "disabled"), (name, st[name])
    assert st["granule"]["status"] == "skipped" and "central_element" in st["granule"]["reason"]
    assert st["coloc"]["status"] == "skipped"
    # the puncta channel resolved by position (laser-line name): a QC flag, not a failure
    assert "channel:foci:name_unmatched_used_index_1" in res["qc_flags"] and res["qc_pass"]

    nuclei = pd.read_csv(out / f"{res['image_id']}__nuclei.csv")
    assert list(nuclei.columns[:len(schema.NUCLEI)]) == schema.NUCLEI
    assert {"dna_mean_intensity", "foci_mean_intensity", "foci_integrated_intensity"} <= set(nuclei.columns)
    assert set(schema.SHARED_META) <= set(nuclei.columns)
    assert len(nuclei) == res["n_nuclei"]
    # no germline stage: every nucleus is the analysis set (a real True, not NA)
    assert nuclei["in_germline"].dtype == bool and nuclei["in_germline"].all()
    # no axis stage: positions are NA on every row
    assert nuclei["axis_position_norm"].isna().all() and nuclei["axis_position_um"].isna().all()
    assert nuclei["n_spots"].notna().all() and nuclei["n_spots"].sum() > 0

    spots = pd.read_csv(out / f"{res['image_id']}__spots.csv")
    assert list(spots.columns) == schema.SPOTS + schema.SHARED_META
    assert len(spots) == int(nuclei["n_spots"].sum()) > 0
    assert (spots["marker"] == "foci").all() and (spots["detector"] == "spotmax").all()
    assert set(spots["nucleus_id"]) <= set(nuclei["nucleus_id"])

    summary = pd.read_csv(out / f"{res['image_id']}__image_summary.csv")
    assert int(summary["n_germline_nuclei"].iloc[0]) == res["n_nuclei"]
    assert summary["mean_spots"].iloc[0] > 0 and summary["total_germline_length_um"].isna().all()
    # a non-worm filename yields the metadata defaults rather than an error
    assert nuclei["image_id"].iloc[0] == "plate1_wellA01_site3"
    assert nuclei["genotype"].iloc[0] == "unknown" and nuclei["sex"].iloc[0] == "unknown"

    # every declared table is on disk (empty for the stages that did not run), plus the images
    for name in schema.TABLES:
        assert (out / f"{res['image_id']}__{name}.csv").exists(), name
    assert pd.read_csv(out / f"{res['image_id']}__granules.csv").empty
    assert pd.read_csv(out / f"{res['image_id']}__coloc.csv").empty
    for suffix in ("montage.png", "nuclei_labels.tif", "spots.tif", "stages.json", "done.json"):
        assert (out / f"{res['image_id']}__{suffix}").exists(), suffix
    man = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["config_hash"] == cfg.hash
    assert not set(WORM_ONLY) & set(man["enabled_stages"])
    assert {"segment", "spots", "write"} <= set(man["enabled_stages"])


def test_germline_off_keeps_every_nucleus_as_a_real_true(tmp_path, patched_reader):  # noqa: F811
    """--no-germline on the worm config: in_germline is True on every row (the whole table is the
    analysis set) and nothing is excluded. The default path (germline on) still writes False rows for
    the off-gonad junk; that is test_end_to_end's and test_golden's territory."""
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("germline.enabled", False)
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["germline"] == {**res["stages"]["germline"], "status": "skipped", "reason": "disabled"}
    nuclei = res["tables"]["nuclei"]
    assert nuclei["in_germline"].dtype == bool and nuclei["in_germline"].all()
    assert len(nuclei) == res["n_nuclei"] and res["n_nuclei"] > 0
    assert int(res["tables"]["image_summary"]["n_germline_nuclei"].iloc[0]) == res["n_nuclei"]
    assert not any(f.startswith("germline:") for f in res["qc_flags"])
    on_disk = pd.read_csv(tmp_path / "r" / f"{res['image_id']}__nuclei.csv")
    assert on_disk["in_germline"].dtype == bool and on_disk["in_germline"].all()
