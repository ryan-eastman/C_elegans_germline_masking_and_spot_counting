"""The signal_qc stage (germquant.signal_qc): per-channel background / noise / SNR / saturation metrics
and the pairwise channel cross-talk table, checked against hand-computed answers on a synthetic stack,
then wired through process_image, the flags switch, the writer and `collect`."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: synthetic 3-channel gonad)

from germquant import batch, pipeline, schema
from germquant.config import load_config
from germquant.signal_qc import (
    REGIONS,
    region_masks,
    run_signal_qc,
    saturation_level,
    signal_flags,
    subsample,
)
from germquant.stages import apply_cli_switches, cli_switches

SPACING = (0.4, 0.2, 0.2)
STRIDE = 2
ROLES = {"dna": 0, "central_element": 1, "foci": 2, "granule": None}   # an unresolved role is skipped


def _stack():
    """(C, Z, Y, X) float32 stack plus labels whose signal metrics can be worked out by hand.

    Channel 0 (dna): a background that alternates 90 / 110 by SUB-SAMPLED column, a cube nucleus at 2000
    and three nucleus voxels at 4095 (12-bit saturation). The cube is centred on the grid, so whatever
    symmetric region the distance transform cuts out of the background holds equal numbers of 90s and
    110s: the background median is exactly 100 and its MAD exactly 10. Channel 1 (central_element) =
    0.5 * channel 0 + 10, a scaled copy (Pearson r = 1, slope = 0.5 in every region). Channel 2 (foci):
    independent Gaussian noise (r near 0), the picture of a channel with no live fluorophore."""
    Z, Y, X = 12, 64, 64
    zz, yy, xx = np.indices((Z, Y, X))
    a = np.where((xx // STRIDE) % 2 == 0, 90.0, 110.0).astype(np.float32)
    labels = np.zeros((Z, Y, X), np.int32)
    labels[4:8, 24:40, 24:40] = 1                 # sub-sampled: z 4..7 of 0..11, rows / cols 12..19 of 0..31
    a[labels > 0] = 2000.0
    for z, y, x in ((5, 30, 30), (5, 30, 32), (6, 28, 26)):   # on the sub-sample grid (even rows / cols)
        a[z, y, x] = 4095.0
    b = (0.5 * a + 10.0).astype(np.float32)
    c = np.random.default_rng(7).normal(500.0, 20.0, a.shape).astype(np.float32)
    return np.stack([a, b, c]), labels


def test_saturation_level_rule():
    assert saturation_level(np.uint16) == 65535.0
    assert saturation_level(np.uint8) == 255.0
    assert saturation_level(np.float32, 12) == 4095.0                       # a known bit depth wins
    assert saturation_level(np.uint16, 12) == 4095.0
    assert saturation_level(np.float32, None, "uint16") == 65535.0           # the file's integer dtype
    assert saturation_level(np.float32) == float(np.finfo(np.float32).max)   # floats: never saturated
    assert saturation_level(np.float32, float("nan")) == float(np.finfo(np.float32).max)
    assert np.isnan(saturation_level("not a dtype"))


def test_subsample_and_regions():
    data, labels = _stack()
    sub = subsample(labels, STRIDE)
    assert sub.shape == (12, 32, 32) and int((sub > 0).sum()) == 4 * 8 * 8
    reg = region_masks(sub, (0.4, 0.8, 0.8), tissue_dilation_um=2.0, ring_um=1.0)
    assert set(reg) == {"nuclei", "perinuclear_ring", "cytoplasm", "tissue", "outside"}
    n = {k: int(v.sum()) for k, v in reg.items()}
    assert n["nuclei"] == 256 and n["perinuclear_ring"] > 0 and n["cytoplasm"] > 0
    assert n["nuclei"] + n["perinuclear_ring"] + n["cytoplasm"] == n["tissue"]
    assert n["tissue"] + n["outside"] == sub.size
    assert not (reg["perinuclear_ring"] & reg["nuclei"]).any() and not (reg["cytoplasm"] & reg["perinuclear_ring"]).any()
    # anisotropy-aware: the ring reaches 1 um = 2.5 planes in z but 1.25 columns at 0.8 um xy spacing
    ring = reg["perinuclear_ring"]
    assert ring[2, 16, 16] and not ring[1, 16, 16]          # 2 planes above the cube (0.8 um) in, 3 (1.2 um) out
    assert ring[6, 16, 11] and not ring[6, 16, 10]          # 1 column (0.8 um) in, 2 columns (1.6 um) out
    empty = region_masks(np.zeros_like(sub), (0.4, 0.8, 0.8), tissue_dilation_um=2.0, ring_um=1.0)
    assert not empty["tissue"].any() and empty["outside"].all()


def test_metrics_match_hand_computation():
    data, labels = _stack()
    acq = {"exp_ms_dna": 10.0, "laser_pow_dna": 50.0, "exp_ms_foci": float("nan"), "laser_pow_foci": 3.0}
    res = run_signal_qc(data, labels, SPACING, ROLES, acq_fields=acq, stride=STRIDE, bit_depth=12,
                        image_id="img")
    assert res["roles"] == ["dna", "central_element", "foci"] and res["n_sampled"] == 12 * 32 * 32
    s = res["summary"]
    assert set(s) == {f"sig_{r}_{m}" for r in res["roles"] for m in schema.SIGNAL_QC_METRICS}
    # channel 0: offset = 0.5th percentile = 90, background median 100, MAD 10 -> noise 14.826,
    # p99 inside the tissue = the nucleus value 2000, SNR = 1900 / 14.826
    assert s["sig_dna_offset"] == 90.0 and s["sig_dna_bg"] == 100.0 and s["sig_dna_bg_source"] == "outside_tissue"
    assert s["sig_dna_noise"] == pytest.approx(14.826, abs=1e-9)
    assert s["sig_dna_p99"] == 2000.0
    assert s["sig_dna_snr"] == pytest.approx(1900.0 / 14.826, rel=1e-9)
    assert s["sig_dna_sat_frac"] == pytest.approx(3 / (12 * 32 * 32), rel=1e-12)      # three voxels at 4095
    assert s["sig_dna_bg_per_ms"] == pytest.approx((100.0 - 90.0) / 10.0) and s["sig_dna_p99_per_ms"] == pytest.approx(191.0)
    assert s["sig_dna_bg_per_ms_pct"] == pytest.approx(1.0 / 50.0)
    # channel 1 is 0.5 * channel 0 + 10: every level scales, the SNR is the same, nothing saturates
    assert s["sig_central_element_offset"] == 55.0 and s["sig_central_element_bg"] == 60.0
    assert s["sig_central_element_p99"] == 1010.0 and s["sig_central_element_noise"] == pytest.approx(7.413, abs=1e-9)
    assert s["sig_central_element_snr"] == pytest.approx(s["sig_dna_snr"], rel=1e-9)
    assert s["sig_central_element_sat_frac"] == 0.0
    # exposure / power unknown for that role -> the normalized levels are NaN (power alone is not enough)
    for k in ("bg_per_ms", "p99_per_ms", "bg_per_ms_pct"):
        assert np.isnan(s[f"sig_central_element_{k}"]) and np.isnan(s[f"sig_foci_{k}"])
    # channel 2 is pure noise: the 99th percentile of a Gaussian sits 2.33 sigma up, so SNR is about 2.3
    assert 1.8 < s["sig_foci_snr"] < 2.9 and abs(s["sig_foci_noise"] - 20.0) < 2.0
    assert "sig_granule_bg" not in s                                                   # unresolved role: no fields
    # the saturation level: the acquisition record's bit depth beats the config key, the dtype rule is last
    assert res["sat_level"] == {"dna": 4095.0, "central_element": 4095.0, "foci": 4095.0}
    res16 = run_signal_qc(data, labels, SPACING, ROLES, acq_fields={"bit_depth_dna": 16}, stride=STRIDE, bit_depth=12)
    assert res16["sat_level"]["dna"] == 65535.0 and res16["sat_level"]["foci"] == 4095.0
    assert res16["summary"]["sig_dna_sat_frac"] == 0.0
    res_f = run_signal_qc(data, labels, SPACING, ROLES, stride=STRIDE)                 # float32 data, no bit depth
    assert res_f["summary"]["sig_dna_sat_frac"] == 0.0 and res_f["sat_level"]["dna"] == float(np.finfo(np.float32).max)
    res_u = run_signal_qc(data, labels, SPACING, ROLES, stride=STRIDE, source_dtype="uint16")
    assert res_u["sat_level"]["dna"] == 65535.0


def test_background_fallback_to_p5():
    data, labels = _stack()
    res = run_signal_qc(data, labels, SPACING, ROLES, stride=STRIDE, min_bg_voxels=10 ** 9)
    s = res["summary"]
    assert s["sig_dna_bg_source"] == "p5" and s["sig_dna_bg"] == 90.0          # 5th percentile of the channel
    assert s["sig_dna_noise"] == 0.0 and np.isnan(s["sig_dna_snr"])            # the darkest 5 percent are all 90
    assert s["sig_dna_offset"] == 90.0 and s["sig_dna_p99"] == 2000.0           # the other metrics are untouched


def test_crosstalk_table_scaled_copy_and_independent_noise():
    data, labels = _stack()
    res = run_signal_qc(data, labels, SPACING, ROLES, stride=STRIDE, image_id="img")
    t = res["crosstalk"]
    assert list(t.columns) == ["image_id", *schema.CHANNEL_CROSSTALK]
    assert (t["image_id"] == "img").all()
    pairs = {(a, b) for a, b in zip(t["role_a"], t["role_b"])}
    assert pairs == {("dna", "central_element"), ("dna", "foci"), ("central_element", "foci")}
    assert len(t) == 3 * len(REGIONS) and list(t["region"].iloc[:3]) == list(REGIONS)
    nuc = t[t["region"] == "nuclei"]
    assert (nuc["n_voxels"] == 256).all()
    copy = t[(t["role_a"] == "dna") & (t["role_b"] == "central_element")].set_index("region")
    for region in REGIONS:
        assert copy.loc[region, "pearson_r"] == pytest.approx(1.0, abs=1e-9), region
        assert copy.loc[region, "slope"] == pytest.approx(0.5, abs=1e-9), region
    noise = t[(t["role_a"] == "dna") & (t["role_b"] == "foci")].set_index("region")
    assert abs(noise.loc["perinuclear_ring", "pearson_r"]) < 0.1 and abs(noise.loc["cytoplasm", "pearson_r"]) < 0.1
    assert abs(noise.loc["nuclei", "pearson_r"]) < 0.25           # 256 voxels, 3 of them spikes
    # a constant channel over a region has no correlation to report
    const = data.copy()
    const[1] = 7.0
    rc = run_signal_qc(const, labels, SPACING, ROLES, stride=STRIDE)["crosstalk"]
    row = rc[(rc["role_a"] == "dna") & (rc["role_b"] == "central_element")]
    assert row["pearson_r"].isna().all() and row["slope"].isna().all() and (row["n_voxels"] > 0).all()


def test_flags_thresholds():
    data, labels = _stack()
    res = run_signal_qc(data, labels, SPACING, ROLES, stride=STRIDE, bit_depth=12)
    s, t, roles = res["summary"], res["crosstalk"], res["roles"]
    # code-default thresholds: the scaled copy is cross-talk in every region; the pure-noise channel
    # (SNR about 2.3) is below min_snr 3; three saturated voxels in 12288 are under max_sat_frac 0.001
    assert signal_flags(s, t, roles) == [
        "signal:foci_low_snr",
        "crosstalk:dna_central_element_nuclei", "crosstalk:dna_central_element_perinuclear_ring",
        "crosstalk:dna_central_element_cytoplasm"]
    fl = signal_flags(s, t, roles, min_snr=200.0, max_sat_frac=1e-4, max_r=2.0)
    assert fl == ["signal:dna_low_snr", "signal:dna_saturated", "signal:central_element_low_snr", "signal:foci_low_snr"]
    assert signal_flags(s, t, roles, min_snr=0.0, max_sat_frac=1.0, max_r=2.0) == []
    # NaN metrics never fire
    nan_summary = {f"sig_dna_{m}": float("nan") for m in schema.SIGNAL_QC_METRICS}
    assert signal_flags(nan_summary, t.assign(pearson_r=np.nan), ["dna"], min_snr=1e9, max_sat_frac=-1.0, max_r=-1.0) == []


def _cfg(**switches):
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    for k, v in switches.items():
        cfg.set(k, v)
    return cfg


def test_stage_runs_by_default_and_adds_only_new_columns(tmp_path, patched_reader):  # noqa: F811
    cfg = _cfg()
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["signal_qc"]["status"] == "ran" and res["stages"]["signal_qc"]["flags"] == []
    assert not any(f.startswith(("signal:", "crosstalk:")) for f in res["qc_flags"])     # flags are off by default
    out = tmp_path / "r"
    isum = pd.read_csv(out / f"{res['image_id']}__image_summary.csv")
    cols = list(isum.columns)
    assert all(c in cols for c in schema.SIGNAL_QC)
    assert cols.index("sig_dna_offset") > cols.index("tail_n_granules_scored")           # appended after the old ones
    row = isum.iloc[0]
    for role in ("dna", "central_element", "foci"):                                         # the three resolved roles
        assert np.isfinite(row[f"sig_{role}_bg"]) and np.isfinite(row[f"sig_{role}_offset"])
        assert row[f"sig_{role}_bg_source"] == "outside_tissue" and row[f"sig_{role}_sat_frac"] == 0.0
        assert np.isnan(row[f"sig_{role}_bg_per_ms"])                                      # acquisition stage off
    for role in ("granule", "lamin"):                                                        # roles the image lacks
        assert all(pd.isna(row[f"sig_{role}_{m}"]) for m in schema.SIGNAL_QC_METRICS)
    assert np.isfinite(row["sig_dna_snr"]) and row["sig_dna_snr"] > 3.0
    # the cross-talk table: 3 pairs x 3 regions, the declared columns first, then the shared block
    ct = pd.read_csv(out / f"{res['image_id']}__channel_crosstalk.csv")
    assert list(ct.columns[:len(schema.CHANNEL_CROSSTALK)]) == schema.CHANNEL_CROSSTALK
    assert set(schema.SHARED_META) <= set(ct.columns) and (ct["image_id"] == res["image_id"]).all()
    assert len(ct) == 9 and set(ct["region"]) == set(REGIONS)
    assert (ct["n_voxels"] > 0).all()
    assert (out / f"{res['image_id']}__channel_crosstalk.parquet").exists()
    assert "channel_crosstalk" in res["tables"] and len(res["tables"]["channel_crosstalk"]) == 9
    # the stage has a provenance sub-hash like the others
    import json

    rec = json.loads((out / f"{res['image_id']}__stages.json").read_text(encoding="utf-8"))
    assert "signal_qc" in rec["stage_hashes"] and "signal_qc" in rec["enabled_stages"]


def test_stage_switch_and_flags_profile(tmp_path, patched_reader):  # noqa: F811
    # switched off: skipped, NaN columns, an empty table still written
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", _cfg(**{"signal_qc.enabled": False}),
                                 tmp_path / "off", prov=None)
    assert res["stages"]["signal_qc"]["status"] == "skipped" and res["stages"]["signal_qc"]["reason"] == "disabled"
    isum = pd.read_csv(tmp_path / "off" / f"{res['image_id']}__image_summary.csv")
    assert isum[schema.SIGNAL_QC].isna().all().all()
    assert pd.read_csv(tmp_path / "off" / f"{res['image_id']}__channel_crosstalk.csv").empty
    # a profile that enables the flags with thresholds every channel trips
    cfg = _cfg(**{"signal_qc.flags.enabled": True, "signal_qc.max_sat_frac": -1.0, "signal_qc.crosstalk.max_r": -1.0,
                  "signal_qc.stride": 2})
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "on", prov=None)
    fl = res["qc_flags"]
    assert {"signal:dna_saturated", "signal:central_element_saturated", "signal:foci_saturated"} <= set(fl)
    assert any(f.startswith("crosstalk:dna_central_element_") for f in fl)
    assert set(res["stages"]["signal_qc"]["flags"]) == {f for f in fl if f.startswith(("signal:", "crosstalk:"))}
    assert res["qc_pass"]                                      # advisory flags, not a hard failure
    # the generated CLI switch
    assert "signal_qc" in [s.name for s in cli_switches()]

    class A:
        no_signal_qc = True

    cfg2 = _cfg()
    h = cfg2.hash
    assert apply_cli_switches(cfg2, A()) == ["skipping the signal_qc stage (--no-signal_qc)"]
    assert cfg2.get("signal_qc.enabled") is False and cfg2.hash != h


def test_collect_stacks_the_crosstalk_table(tmp_path, patched_reader):  # noqa: F811
    cfg = _cfg()
    root = tmp_path / "res"
    for iid in ("20251105_n2_nohs_HERM_001.nd2", "20251105_n2_hs_MALE_002.nd2"):
        pipeline.process_image(iid, cfg, root / iid[:-4], prov=None)
    counts = batch.collect(root)
    assert counts["channel_crosstalk"] == 18
    stacked = pd.read_csv(root / "batch_channel_crosstalk.csv")
    assert set(stacked["image_id"]) == {"20251105_n2_nohs_HERM_001", "20251105_n2_hs_MALE_002"}
    assert list(stacked.columns[:len(schema.CHANNEL_CROSSTALK)]) == schema.CHANNEL_CROSSTALK
