"""The live GFP::SYP-2 profile (config/profiles/live_syp2.yaml with config/channel_maps/live_gfp_syp2.yaml) and
the segmentation.nuclei.channel key it relies on (docs/TIMELAPSE.md, "Segmenting without DAPI").

The key is live: the segment stage segments the channel whose ROLE it names, with the code default "dna", so
every frozen recipe resolves the pixels it always did (test_golden gates that). The live map has no dna role at
all: central_element (EGFP) is segmented, brightfield is a new role no stage consumes, and the stages that need
RAD-51, PGL-1 or lamin are off or gate themselves off. A synthetic two-channel series (nuclear GFP blobs with an
SC bundle, a transmitted-light gradient) with a planted drift runs through process_series end to end on the CPU.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel fixed gonad)
from test_segment_backends import registry  # noqa: F401  (fixture: a private backend registry)

from germquant import cli, pipeline, schema
from germquant.config import load_config
from germquant.io.channel_map import ChannelMap
from germquant.io.nd2_reader import Stack
from germquant.segment import backends as B
from germquant.stages import stage_enabled
from germquant.timelapse import process_series

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "config" / "profiles" / "live_syp2.yaml"
CHANNEL_MAP = ROOT / "config" / "channel_maps" / "live_gfp_syp2.yaml"
SPACING = (0.4, 0.2, 0.2)
LIVE_CHANNELS = ["Brightfield", "EGFP"]          # the October 2026 series, in its channel order
ND2 = "20261005_DLW93_nohs_gfpsyp2_5min_2hr_herm.nd2"
STEM = "20261005_DLW93_nohs_gfpsyp2_5min_2hr_herm"
T = 4
PERIOD_S = 300.0
DRIFT_PX = (1, 3)                                # (y, x) pixels of drift per timepoint
OFF = ("spots", "granule", "coloc", "envelope", "partition", "staging")
ON = ("segment", "measure", "signal_qc", "germline", "axis", "qc", "render", "write")


def _live_stack() -> Stack:
    """A two-channel live stack: 14 GFP::SYP-2 nuclei along x (a nuclear Gaussian plus a bright SC bundle
    through each), a transmitted-light channel that is a flat field with a gradient, and a little noise so
    the signal QC has a background sigma to measure. No DNA channel."""
    Z, Y, X = 20, 80, 320
    rng = np.random.default_rng(7)
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    gfp = np.zeros((Z, Y, X), np.float32)
    for cx in range(20, 300, 20):
        r2 = (((zz - 10) * SPACING[0]) ** 2 + ((yy - 40) * SPACING[1]) ** 2 + ((xx - cx) * SPACING[2]) ** 2)
        gfp += np.exp(-r2 / (2 * 1.2 ** 2)) * 1500.0
        gfp[10, 40, max(0, cx - 5):min(X, cx + 5)] += 1500.0
    gfp += rng.normal(100.0, 5.0, gfp.shape).astype(np.float32)
    bf = (2000.0 + 3.0 * xx + rng.normal(0.0, 20.0, gfp.shape)).astype(np.float32)
    return Stack(data=np.stack([bf, gfp]), spacing=SPACING, channel_names=LIVE_CHANNELS, path="x.nd2", spacing_ok=True)


def live_series_stacks(n_t: int = T) -> list[Stack]:
    base = _live_stack()
    _, _, Y, X = base.data.shape
    out = []
    for k in range(n_t):
        dy, dx = DRIFT_PX[0] * k, DRIFT_PX[1] * k
        data = np.zeros_like(base.data)
        data[:, :, dy:, dx:] = base.data[:, :, : Y - dy, : X - dx]        # whole-pixel shift, zero fill
        out.append(Stack(data=data, spacing=SPACING, channel_names=LIVE_CHANNELS, path="x.nd2", spacing_ok=True))
    return out


def live_series_meta(n_t: int = T) -> dict:
    return {"sizes": {"T": n_t, "C": 2, "Z": 20, "Y": 80, "X": 320}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": LIVE_CHANNELS, "dtype": "uint16", "is_2d": False,
            "n_timepoints": n_t, "period_s": PERIOD_S, "t_seconds": [k * PERIOD_S for k in range(n_t)]}


def patch_live_series(monkeypatch, n_t: int = T) -> list:
    """Point the pipeline's reader seam at a synthetic live series; returns the list of t values read."""
    stacks = live_series_stacks(n_t)
    meta = live_series_meta(n_t)
    calls: list = []

    def read_stack(path, *, xy_stride=1, z_range=None, as_float=False, t=None):
        calls.append(t)
        k = 0 if t is None else int(t)
        if not 0 <= k < n_t:
            raise ValueError(f"t={k} out of range")
        return stacks[k]

    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", read_stack)
    return calls


@pytest.fixture
def patched_live_series(monkeypatch):
    return patch_live_series(monkeypatch, T)


def live_cfg():
    cfg = load_config(PROFILE)
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"      # CPU, deterministic, as every e2e test
    return cfg


# ---- the channel map and the profile ------------------------------------------------------------------

def test_live_channel_map_has_no_dna_role_and_resolves_the_series_channels():
    cm = ChannelMap.from_yaml(CHANNEL_MAP)
    assert list(cm.roles) == ["central_element", "brightfield"]
    assert not cm.has_role("dna") and cm.roles["central_element"].required and not cm.roles["brightfield"].required
    assert cm.marker("central_element") == "GFP::SYP-2" and cm.marker("brightfield") == "brightfield"
    assert cm.resolve(LIVE_CHANNELS) == ({"central_element": 1, "brightfield": 0}, [])
    assert cm.resolve(["TD", "GFP"]) == ({"central_element": 1, "brightfield": 0}, [])
    # a GFP-only file: brightfield never resolves by position (no collision), no flag
    assert cm.resolve(["488 nm"]) == ({"central_element": 0, "brightfield": None}, [])
    # names that match nothing: central_element falls back to channel 1 and says so; a DAPI file has no SC channel
    roles, flags = cm.resolve(["ch0", "ch1"])
    assert roles["central_element"] == 1 and "channel:central_element:name_unmatched_used_index_1" in flags
    roles, flags = cm.resolve(["DAPI"])
    assert roles["central_element"] is None and "channel:central_element:MISSING_REQUIRED" in flags


def test_live_profile_switches_and_the_base_hash_is_untouched():
    sentinel = (ROOT / "tests" / "golden" / "synthetic" / "config_hash.txt").read_text().strip()
    base = load_config(ROOT / "config" / "config.yaml")
    assert base.hash == sentinel
    assert base.get("segmentation.nuclei.channel") == "dna"                 # the frozen recipes name dna
    prof = load_config(PROFILE)
    assert prof.hash != base.hash and prof.hash == load_config(PROFILE).hash
    assert load_config(ROOT / "config" / "config.yaml").hash == sentinel
    assert prof.get("io.channel_map") == "config/channel_maps/live_gfp_syp2.yaml"
    assert set(prof.channel_map.roles) == {"central_element", "brightfield"} and not prof.channel_map.has_role("dna")
    assert prof.get("segmentation.nuclei.channel") == "central_element"
    assert prof.get("segmentation.nuclei.cellpose_model") == base.get("segmentation.nuclei.cellpose_model")  # inherited
    assert prof.get("signal_qc.enabled") is True and prof.get("signal_qc.flags.enabled") is False
    assert prof.get("timelapse.keep_reference_volume") is True
    assert base.get("timelapse.keep_reference_volume") is None             # a new key, read with a code default
    for name in OFF:
        assert stage_enabled(prof, name) is False, name
    for name in ON:
        assert stage_enabled(prof, name) is True, name
    assert stage_enabled(prof, "audit") is False and stage_enabled(prof, "granule_tail") is False


# ---- the channel key ----------------------------------------------------------------------------------

def _spy_backend(seen: dict):
    def spy(dna, spacing, *, model, diameter_um, device, cfg_section):
        from germquant.segment import nuclei as N

        seen["array"] = np.asarray(dna).copy()
        seen["model"] = model
        return N._classical(dna, spacing, diameter_um)
    return spy


def test_segmentation_channel_key_picks_the_role(tmp_path, patched_reader, registry):  # noqa: F811
    """segmentation.nuclei.channel names the ROLE the segment stage segments: dna by default (the goldens),
    any other resolved role on request, and an unresolved role falls back to channel 0 with a flag that
    names the role (the dna text is the flag the fixed path always wrote)."""
    seen: dict = {}
    B.register_backend("spy", _spy_backend(seen))
    stack = pipeline.read_stack("x")                                       # the fixture's synthetic stack

    def run(channel, sub):
        cfg = load_config(ROOT / "config" / "config.yaml")
        cfg._data["segmentation"]["nuclei"]["method"] = "auto"
        cfg._data["segmentation"]["nuclei"]["backend"] = "spy"
        if channel is not None:
            cfg._data["segmentation"]["nuclei"]["channel"] = channel
        return pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / sub, prov=None)

    res = run(None, "dna_default")                                         # config.yaml says dna
    assert np.array_equal(seen["array"], stack.data[0]) and res["stages"]["segment"]["status"] == "ran"
    assert not any(f.startswith("segment:no_") for f in res["qc_flags"])
    res = run("foci", "foci")                                              # the RAD-51 channel is segmented
    assert np.array_equal(seen["array"], stack.data[2])
    assert not np.array_equal(seen["array"], stack.data[0])
    res = run("granule", "granule")                                        # unresolved in a 3-channel image
    assert np.array_equal(seen["array"], stack.data[0])
    assert "segment:no_granule_channel_using_first_channel" in res["qc_flags"]
    isum = res["tables"]["image_summary"]
    assert isum["segmentation_method"].iloc[0] == "spy"


# ---- the synthetic live series end to end -------------------------------------------------------------

def test_live_series_runs_end_to_end(tmp_path, patched_live_series):
    cfg = live_cfg()
    out = tmp_path / "live"
    res = process_series(ND2, cfg, out, timepoints="all")
    assert res["done"] == list(range(T)) and res["failed"] == [] and patched_live_series == list(range(T))

    centroids = {}
    for k in range(T):
        sub, iid = out / f"t{k:03d}", f"{STEM}_t{k:03d}"
        rec = json.loads((sub / f"{iid}__stages.json").read_text(encoding="utf-8"))["stages"]
        assert not [s for s, o in rec.items() if o["status"] == "failed"], rec
        for name in ON:
            assert rec[name]["status"] == "ran", (name, rec[name])
        for name in OFF:
            assert (rec[name]["status"], rec[name]["reason"]) == ("skipped", "disabled"), (name, rec[name])
        nuclei = pd.read_csv(sub / f"{iid}__nuclei.csv")
        assert len(nuclei) >= 10 and nuclei["in_germline"].astype(bool).all()        # one connected GFP-positive row
        # intensity columns for the roles present; the dna ones are there and NaN (no DNA stain)
        for role in ("central_element", "brightfield"):
            assert nuclei[f"{role}_mean_intensity"].notna().all(), role
        for m in ("mean", "max", "integrated"):
            assert f"dna_{m}_intensity" in nuclei.columns and nuclei[f"dna_{m}_intensity"].isna().all()
        assert nuclei["central_element_mean_intensity"].median() > 500 > nuclei["brightfield_mean_intensity"].std()
        assert not nuclei["dapi_present"].astype(bool).any()
        assert json.loads(nuclei["channel_map_json"].iloc[0]) == {"central_element": 1, "brightfield": 0}
        assert nuclei["genotype"].iloc[0] == "DLW93" and nuclei["sex"].iloc[0] == "herm" and nuclei["treatment"].iloc[0] == "control"
        assert nuclei["axis_position_um"].notna().any()
        isum = pd.read_csv(sub / f"{iid}__image_summary.csv")
        assert isum["segmentation_method"].iloc[0] == "classical"
        assert isum["sig_central_element_snr"].iloc[0] > 3 and isum["sig_brightfield_bg"].notna().all()
        assert "sig_dna_snr" in isum.columns and isum["sig_dna_snr"].isna().all()      # declared, no such channel
        assert pd.read_csv(sub / f"{iid}__spots.csv").empty and pd.read_csv(sub / f"{iid}__granules.csv").empty
        assert "spots:disabled_segmentation_only" in str(isum["qc_flags"].iloc[0])
        assert (sub / f"{iid}__montage.png").exists() and (sub / f"{iid}__nuclei_labels.tif").exists()
        man = json.loads((sub / "run_manifest.json").read_text(encoding="utf-8"))
        assert man["config"]["segmentation"]["nuclei"]["channel"] == "central_element"
        assert man["config"]["timelapse"]["keep_reference_volume"] is True and man["t_index"] == k
        centroids[k] = nuclei.loc[nuclei["in_germline"].astype(bool), ["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(float)

    # the planted drift is in the numbers: nearest-neighbour matched nucleus centroids move by k * (1, 3) pixels
    from scipy.spatial import cKDTree

    ref = centroids[0]
    for k in range(1, T):
        pts = centroids[k]
        _, idx = cKDTree(ref).query(pts)
        step = np.median(pts - ref[idx], axis=0)
        assert step[2] == pytest.approx(k * DRIFT_PX[1] * SPACING[2], abs=0.1), ("x", k, step)
        assert step[1] == pytest.approx(k * DRIFT_PX[0] * SPACING[1], abs=0.1), ("y", k, step)
        assert step[0] == pytest.approx(0.0, abs=0.1), ("z", k, step)

    # the series summary carries the signal QC of both live roles, in channel-map order, and nothing for dna
    ss = pd.read_csv(out / f"{STEM}__series_summary.csv")
    assert list(ss["t_index"]) == list(range(T)) and ss["n_nuclei"].ge(10).all() and ss["qc_pass"].all()
    assert [c for c in ss.columns if c.startswith("sig_")] == ["sig_central_element_bg", "sig_central_element_snr",
                                                                "sig_brightfield_bg", "sig_brightfield_snr"]
    assert list(ss.columns[:len(schema.SERIES_SUMMARY)]) == schema.SERIES_SUMMARY
    assert (ss["genotype"] == "DLW93").all() and (ss["config_hash"] == cfg.hash).all()


def test_live_profile_on_the_command_line(tmp_path, patched_live_series):
    """`--profile live_syp2` is found by name; the run itself goes through a study config that extends the
    profile with the classical segmenter (the shipped profile says `method: auto`, which would run Cellpose-SAM
    on the CPU here, as the fixed-image CLI tests avoid in the same way)."""
    from germquant.config import find_profile

    assert find_profile("live_syp2") == PROFILE
    study = tmp_path / "live_classical.yaml"
    study.write_text(f"extends: {PROFILE.as_posix()}\nsegmentation:\n  nuclei:\n    method: classical\n", encoding="utf-8")
    out = tmp_path / "cli"
    assert cli.main(["run", ND2, "--config", str(study), "--out", str(out), "--t", "0-1"]) == 0
    assert sorted(p.name for p in out.iterdir() if p.is_dir()) == ["t000", "t001"]
    assert patched_live_series == [0, 1]
    ss = pd.read_csv(out / f"{STEM}__series_summary.csv")
    assert len(ss) == 2 and ss["qc_pass"].all()
    assert [c for c in ss.columns if c.startswith("sig_")] == ["sig_central_element_bg", "sig_central_element_snr",
                                                                "sig_brightfield_bg", "sig_brightfield_snr"]
    man = json.loads((out / "t000" / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["config"]["segmentation"]["nuclei"]["channel"] == "central_element"
    assert man["config"]["io"]["channel_map"].endswith("live_gfp_syp2.yaml")
