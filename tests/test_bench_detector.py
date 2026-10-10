"""The detector bench harness (germquant.benchmark.detector, scripts/bench_cache.py,
scripts/bench_detector.py) on synthetic data:

- the chunked .npy writer and the data sha256;
- a cache round trip from a real pipeline run on a synthetic uint16 stack: the cached label image is the
  one the spots stage received, and run_backend calls detect_spots with exactly the pipeline's arguments
  and reproduces its spots table (log backend always, SpotMAX when it is installed);
- the merge rule of params over the cached spots: section is the pipeline's spot instance rule;
- the test guard (library and script), the trial log, the final-run log and its limit;
- the runner end to end on tiny planted-foci stacks with the log backend, scored against truth files.
Nothing here needs microscopy data or a GPU.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import CHANNELS, SPACING, _synthetic_stack

from germquant import pipeline
from germquant.benchmark import detector as bd
from germquant.benchmark.manifest import Entry
from germquant.config import load_config
from germquant.io.nd2_reader import Stack

ROOT = Path(__file__).resolve().parents[1]
SP = (0.2, 0.108, 0.108)        # dz, dy, dx: the real anisotropy
SHAPE = (36, 128, 128)


def _load_script(name: str):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------- the writer
def test_write_npy_chunked_and_sha(tmp_path):
    rng = np.random.default_rng(0)
    arr = rng.integers(0, 60000, size=(8, 5, 7)).astype(np.uint16)
    sha = bd.write_npy_chunked(tmp_path / "a.npy", arr, z_chunk=3)
    back = np.load(tmp_path / "a.npy")
    assert back.dtype == np.uint16 and np.array_equal(back, arr)
    assert sha == hashlib.sha256(arr.tobytes()).hexdigest() == bd.array_sha256(np.load(tmp_path / "a.npy", mmap_mode="r"))
    assert not list(tmp_path.glob("*.partial*"))               # the temporary name is gone
    counts = bd.label_voxel_counts(np.array([[[0, 1, 1], [3, 3, 3]]], np.int32))
    assert counts.tolist() == [1, 2, 0, 3]


def test_effective_section_is_the_pipeline_merge_rule():
    class Shim:
        def __init__(self, spots):
            self.spots = spots

        def get(self, key, default=None):
            return self.spots if key == "spots" else default

    golden = {"effect_size_min": 3.0, "log": {"diameter_um": 0.5, "min_snr": 2.0}, "enabled": True}
    for params in ({}, {"effect_size_min": 2.5}, {"log": {"min_snr": 3.0}}, {"log": {"x": 1}, "new_key": [1, 2]}):
        assert bd.effective_section(golden, params) == pipeline.spot_backend_section(Shim(golden), params)
    assert golden == {"effect_size_min": 3.0, "log": {"diameter_um": 0.5, "min_snr": 2.0}, "enabled": True}
    with pytest.raises(ValueError, match="backend"):
        bd.effective_section(golden, {"backend": "log"})
    assert bd.params_sha256({}) == hashlib.sha256(b"{}").hexdigest()
    assert bd.params_sha256({"b": 1, "a": 2}) == bd.params_sha256({"a": 2, "b": 1})


# ------------------------------------------------------- a real pipeline run as the golden run
@pytest.fixture
def uint16_stack(monkeypatch):
    """The synthetic 3-channel gonad of test_end_to_end as uint16 (the real .nd2 dtype), fed to the
    pipeline through its reader seams."""
    s = _synthetic_stack()
    stack = Stack(data=np.round(s.data).astype(np.uint16), spacing=s.spacing, channel_names=list(CHANNELS),
                  path=s.path, spacing_ok=True)
    meta = {"sizes": {"C": 3, "Z": stack.data.shape[1], "Y": stack.data.shape[2], "X": stack.data.shape[3]},
            "spacing": SPACING, "spacing_ok": True, "channel_names": list(CHANNELS), "dtype": "uint16",
            "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)
    return stack


def _split_record(tag: str, nd2: Path) -> dict:
    data = nd2.read_bytes()
    return {"gonads": [{"tag": tag, "nd2": {"path": str(nd2), "size_bytes": len(data),
                                             "sha256": hashlib.sha256(data).hexdigest()}}]}


@pytest.mark.parametrize("backend", ["log", "spotmax"])
def test_cache_round_trip_reproduces_the_spots_stage(tmp_path, monkeypatch, uint16_stack, backend):
    if backend == "spotmax":
        pytest.importorskip("spotmax")
    import germquant.spots as gs

    calls = []
    real = gs.detect_spots

    def spy(*a, **k):
        calls.append((a, k))
        return real(*a, **k)

    monkeypatch.setattr(gs, "detect_spots", spy)
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    if backend != "spotmax":
        cfg.set("spots.backend", backend)
    nd2 = tmp_path / "20251105_n2_nohs_HERM_001.nd2"
    nd2.write_bytes(b"a synthetic stand-in for the registered stack")
    res = pipeline.process_image(nd2, cfg, tmp_path / "golden", prov=None)
    assert res["stages"]["spots"]["status"] == "ran", res["stages"]["spots"]
    assert len(calls) == 1 and len(res["tables"]["spots"]) > 0

    tag = "synth_golden"
    entry = Entry(tag=tag, results_dir=tmp_path / "golden", split="train", nd2=nd2, stem=res["image_id"],
                  extra={"role": "train", "condition": "noHS", "sex": "herm", "spot_truth_usable": True})
    st = bd.cache_gonad(entry, tmp_path / "cache", split_record=_split_record(tag, nd2),
                        reader=lambda p: uint16_stack, z_chunk=3, verify_sha256=True)
    assert st["status"] == "written", st
    assert bd.cache_gonad(entry, tmp_path / "cache", reader=lambda p: uint16_stack)["status"] == "exists"

    cache = bd.load_cache(tag, tmp_path / "cache", verify=True)
    m = cache.meta
    assert np.array_equal(cache.foci, uint16_stack.data[2]) and cache.foci.dtype == np.uint16
    assert np.array_equal(cache.dna, uint16_stack.data[0]) and cache.dna.dtype == np.uint16
    # the cached labels are the label image the spots stage received (all nuclei), as int32
    assert cache.labels.dtype == np.int32 and np.array_equal(cache.labels, calls[0][0][1])
    assert m["spacing_zyx_um"] == list(SPACING) and cache.spacing == tuple(SPACING)
    assert m["role_indices"] == {"dna": 0, "central_element": 1, "foci": 2, "granule": None}
    golden_manifest = json.loads((tmp_path / "golden" / "run_manifest.json").read_text(encoding="utf-8"))
    assert m["spots_config"] == golden_manifest["config"]["spots"]
    assert m["nd2"]["sha256_verified"] and m["nd2"]["sha256"] == _split_record(tag, nd2)["gonads"][0]["nd2"]["sha256"]
    assert m["marker"] == "RAD-51" and m["role"] == "train"
    nuc = res["tables"]["nuclei"]
    assert m["germline_label_ids"] == sorted(int(i) for i in nuc.loc[nuc["in_germline"].astype(bool), "nucleus_id"])
    assert m["n_nuclei"] == len(nuc) and int(cache.labels.max()) == m["labels"]["max_label"]

    (a0, k0), = calls                                    # what the pipeline handed detect_spots
    calls.clear()
    per_spot, per_nuc = bd.run_backend(tag, backend, {}, cache_root=tmp_path / "cache")
    (a1, k1), = calls                                    # what run_backend handed it
    assert np.array_equal(a1[0], a0[0]) and a1[0].dtype == a0[0].dtype == np.uint16
    assert np.array_equal(a1[1], a0[1]) and a1[1].dtype == a0[1].dtype
    assert tuple(a1[2]) == tuple(a0[2]) == tuple(SPACING)
    assert type(a1[0]) is np.ndarray and type(a1[1]) is np.ndarray      # in memory, not file maps
    assert k1 == k0
    spots_section = dict(cfg.get("spots"))
    assert k1 == {**bd.pipeline_spot_kwargs(spots_section), "marker": "RAD-51", "backend": backend,
                  "cfg_section": spots_section}
    pd.testing.assert_frame_equal(per_spot, res["tables"]["spots"])
    golden = bd.read_golden_spots(tmp_path / "golden", res["image_id"])
    cmp = bd.compare_with_golden(per_spot, per_nuc, golden, pd.read_csv(tmp_path / "golden" / f"{res['image_id']}__nuclei.csv"))
    assert cmp["identical"], cmp
    bad = per_spot.copy()
    bad.loc[0, "x_um"] += 1e-9
    cmp2 = bd.compare_with_golden(bad, per_nuc, golden)
    assert not cmp2["identical"] and cmp2["columns"]["x_um"]["n_different"] == 1


def test_pipeline_kwargs_follow_the_spots_stage(tmp_path, monkeypatch, uint16_stack):
    """The keyword arguments of the spots stage are the ones pipeline_spot_kwargs builds from the same
    section, also for non-default values (log backend, so SpotMAX is not needed): a drift in either
    place fails here."""
    import germquant.spots as gs
    from germquant.spots.detect import detect_spots as real

    calls = []

    def spy(*a, **k):
        calls.append((a, k))
        return real(*a, **k)

    monkeypatch.setattr(gs, "detect_spots", spy)
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("spots.backend", "log")
    cfg.set("spots.effect_size_min", 2.5)
    cfg.set("spots.z_merge_gap_um", 0.6)
    cfg.set("spots.log", {"min_snr": 2.5})
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["spots"]["status"] == "ran"
    (_, k_pipe), = calls
    section = dict(cfg.get("spots"))
    assert {k: v for k, v in k_pipe.items() if k not in ("marker", "backend", "cfg_section")} ==         bd.pipeline_spot_kwargs(section)
    assert k_pipe["cfg_section"] == bd.effective_section(section, {}) and k_pipe["backend"] == "log"
    assert bd.pipeline_spot_kwargs(section)["effect_size_min"] == 2.5


# --------------------------------------------------------------------- planted-foci gonads
def _planted(counts=(2, 4, 6, 3), seed=11, sigma_um=0.17):
    """Four spherical nuclei (radius 1.5 um) with ``counts`` Gaussian foci each, Gaussian noise, uint16."""
    rng = np.random.default_rng(seed)
    sp = np.array(SP)
    zz, yy, xx = np.indices(SHAPE).astype(np.float32)
    zum, yum, xum = zz * sp[0], yy * sp[1], xx * sp[2]
    labels = np.zeros(SHAPE, np.int32)
    img = np.full(SHAPE, 30.0, np.float32)
    centers = [(3.6, 3.6, 3.6), (3.6, 3.6, 10.2), (3.6, 10.2, 3.6), (3.6, 10.2, 10.2)]
    planted = []
    for k, (c, n) in enumerate(zip(centers, counts), start=1):
        c = np.array(c)
        inside = ((zum - c[0]) ** 2 + (yum - c[1]) ** 2 + (xum - c[2]) ** 2) <= 1.5 ** 2
        labels[inside] = k
        img[inside] = 120.0
        pts: list[np.ndarray] = []
        while len(pts) < n:
            v = rng.uniform(-0.9, 0.9, size=3)
            if np.linalg.norm(v) > 0.9 or any(np.linalg.norm(c + v - q) < 0.7 for q in pts):
                continue
            pts.append(c + v)
        for p in pts:
            img += rng.uniform(400, 800) * np.exp(-((zum - p[0]) ** 2 + (yum - p[1]) ** 2 + (xum - p[2]) ** 2)
                                                  / (2 * sigma_um ** 2))
            planted.append((k, *p))
    img += rng.normal(0, 10, size=SHAPE).astype(np.float32)
    img = np.clip(np.round(img), 0, 65535).astype(np.uint16)
    dna = np.where(labels > 0, 1000, 50).astype(np.uint16)
    return img, dna, labels, pd.DataFrame(planted, columns=["nucleus_id", "z_um", "y_um", "x_um"])


def _write_gonad(root: Path, tag: str, role: str, condition: str, counts, seed) -> dict:
    """Golden outputs, the dummy registered stack, the truth files and the cache for one planted gonad."""
    import tifffile

    foci, dna, labels, planted = _planted(counts, seed)
    stem = f"img_{tag}"
    gdir = root / "golden" / tag
    gdir.mkdir(parents=True)
    ids = np.arange(1, labels.max() + 1)
    zz, yy, xx = np.nonzero(labels)
    lab = labels[zz, yy, xx]
    cz, cy, cx = (np.bincount(lab, weights=v)[1:] / np.bincount(lab)[1:] for v in (zz, yy, xx))
    nuclei = pd.DataFrame({"nucleus_id": ids, "centroid_z_um": cz * SP[0], "centroid_y_um": cy * SP[1],
                           "centroid_x_um": cx * SP[2], "in_germline": [True, True, True, False],
                           "n_spots": 0, "voxel_dz_um": SP[0], "voxel_dy_um": SP[1], "voxel_dx_um": SP[2]})
    nuclei.to_csv(gdir / f"{stem}__nuclei.csv", index=False)
    tifffile.imwrite(gdir / f"{stem}__nuclei_labels.tif", labels)
    pd.DataFrame({"spot_id": [0], "nucleus_id": [1], "marker": ["RAD-51"], "z_um": [0.0], "y_um": [0.0],
                  "x_um": [0.0], "effect_size": [5.0], "spot_mean_intensity": [1.0],
                  "detector": ["spotmax"]}).to_csv(gdir / f"{stem}__spots.csv", index=False)
    pd.DataFrame({"image_id": [stem], "channel_map_json": [json.dumps({"dna": 0, "central_element": 1,
                                                                        "foci": 2, "granule": None})]}
                 ).to_csv(gdir / f"{stem}__image_summary.csv", index=False)
    spots_cfg = {"enabled": True, "spot_radius_um": 0.3, "thresholding_method": "threshold_triangle",
                 "effect_size_min": 3.0, "merge_z_columns": True, "z_merge_gap_um": 0.8,
                 "z_merge_valley_frac": 0.8, "max_spot_candidates": 30000}
    (gdir / "run_manifest.json").write_text(json.dumps({"git_sha": "synthetic", "config_hash": "synthetic",
                                                        "config": {"spots": spots_cfg}}), encoding="utf-8")
    nd2 = root / "stacks" / f"{stem}.nd2"
    nd2.parent.mkdir(parents=True, exist_ok=True)
    nd2.write_bytes(f"stand-in for {tag}".encode())

    # truth: the planted foci at their voxels; the first focus of every nucleus is Class B
    tdir = root / "truth" / tag
    tdir.mkdir(parents=True)
    zi = np.rint(planted["z_um"] / SP[0]).astype(int)
    yi = np.rint(planted["y_um"] / SP[1]).astype(int)
    xi = np.rint(planted["x_um"] / SP[2]).astype(int)
    first = ~planted["nucleus_id"].duplicated()
    pd.DataFrame({"spot_id": np.arange(len(planted)) + 100, "x_um": xi * SP[2], "y_um": yi * SP[1],
                  "z_um": zi * SP[0], "xi": xi, "yi": yi, "zi": zi, "radius_um": 0.2, "quality": np.nan,
                  "spot_class": np.where(first, "Class B", "Class A")}).to_csv(tdir / f"{stem}__imaris_spots.csv", index=False)
    pd.DataFrame({"surface_id": ids, "label": ids, "x_um": cx * SP[2], "y_um": cy * SP[1], "z_um": cz * SP[0],
                  "xi": np.rint(cx).astype(int), "yi": np.rint(cy).astype(int), "zi": np.rint(cz).astype(int),
                  "volume_um3": np.bincount(lab)[1:] * np.prod(SP), "n_voxels_mask": np.bincount(lab)[1:]}
                 ).to_csv(tdir / f"{stem}__imaris_nuclei.csv", index=False)
    tifffile.imwrite(tdir / f"{stem}__imaris_labels.tif", labels)
    (tdir / f"{stem}__imaris_params.json").write_text(json.dumps(
        {"voxel_zyx_um_ours": list(SP), "voxel_xyz_um_imaris": [0.108, 0.108, 0.2], "region_growing": False}),
        encoding="utf-8")

    entry = {"tag": tag, "role": role, "split": role, "condition": condition, "sex": "herm", "stem": stem,
             "nd2": str(nd2), "ims": str(root / f"{stem}.ims"), "results_dir": str(gdir), "truth_dir": str(tdir),
             "spot_truth_usable": True}
    stack = Stack(data=np.stack([dna, np.zeros_like(dna), foci]), spacing=SP, channel_names=list(CHANNELS),
                  path=nd2, spacing_ok=True)
    e = Entry(tag=tag, results_dir=gdir, split=role, truth_dir=tdir, nd2=nd2, stem=stem,
              extra={"role": role, "condition": condition, "sex": "herm", "spot_truth_usable": True})
    st = bd.cache_gonad(e, root / "cache", reader=lambda p: stack, z_chunk=5)
    assert st["status"] == "written", st
    return {"entry": entry, "planted": planted}


@pytest.fixture
def planted_bench(tmp_path):
    root = tmp_path
    gonads = [_write_gonad(root, "synth_noHS_a", "train", "noHS", (2, 4, 6, 3), 11),
              _write_gonad(root, "synth_HS_b", "train", "HS", (5, 1, 3, 7), 12),
              _write_gonad(root, "synth_val_c", "val", "noHS", (3, 3, 5, 2), 13),
              _write_gonad(root, "synth_test_d", "test", "HS", (4, 2, 2, 6), 14)]
    unusable = dict(gonads[0]["entry"], tag="synth_rg_e", spot_truth_usable=False)
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({"gonads": [g["entry"] for g in gonads] + [unusable]}), encoding="utf-8")
    return {"root": root, "manifest": manifest, "gonads": gonads}


def _bench_args(b, *extra):
    r = b["root"]
    return ["--manifest", str(b["manifest"]), "--cache-root", str(r / "cache"), "--trials-log",
            str(r / "trials.jsonl"), "--final-log", str(r / "final_runs.jsonl"), "--n-boot", "50", *extra]


def test_runner_log_backend_on_planted_stacks(planted_bench, capsys):
    b = planted_bench
    mod = _load_script("bench_detector.py")
    out = b["root"] / "out_train"
    rc = mod.main(["--backend", "log", "--params", "{}", "--split", "train", "--out", str(out), *_bench_args(b)])
    text = capsys.readouterr().out
    assert rc == 0, text
    per_gonad = pd.read_csv(out / "per_gonad.csv")
    assert per_gonad["tag"].tolist() == ["synth_noHS_a", "synth_HS_b"]          # the unusable gonad is not in train
    assert (per_gonad["role"] == "train").all() and per_gonad["split"].tolist() == ["noHS", "HS"]
    for _, r in per_gonad.iterrows():
        assert r["recall"] >= 0.8 and r["precision"] >= 0.8, r[["tag", "precision", "recall"]].to_dict()
        assert r["traced_rule"] == "overlap" and r["n_traced"] == 3 and r["n_germline_nuclei"] == 3
        assert r["ccc_all"] > 0.8 and r["n_spots_reconcile"]
        assert r["n_expert_spots_class_a"] == r["n_expert_spots"] - 4
    per_nuc = pd.read_csv(out / "per_nucleus.csv")
    assert {"expert_n", "expert_n_class_a", "ours_n", "traced", "role"} <= set(per_nuc.columns)
    assert (per_nuc["expert_n"] - per_nuc["expert_n_class_a"] == 1).all()     # one Class B focus per nucleus
    pooled = json.loads((out / "pooled.json").read_text(encoding="utf-8"))
    for k in ("per_nucleus_all", "per_nucleus_traced", "per_nucleus_all_class_a", "per_nucleus_traced_class_a",
              "spots_all", "spots_without_class_b", "spots_detector_only", "by_split", "skipped"):
        assert k in pooled, k
    assert pooled["n_gonads"] == 2 and pooled["backend"] == "log" and pooled["bootstrap"]["seed"] == bd.BOOT_SEED
    assert set(pooled["by_split"]) == {"noHS", "HS"} and "per_nucleus_all_class_a" in pooled["by_split"]["HS"]
    assert (out / "detections" / "synth_noHS_a__spots.csv").is_file()
    # the val split runs the val gonad only; the region-growing gonad is skipped from train when named
    rc = mod.main(["--backend", "log", "--split", "val", "--out", str(b["root"] / "out_val"), *_bench_args(b)])
    assert rc == 0
    assert pd.read_csv(b["root"] / "out_val" / "per_gonad.csv")["tag"].tolist() == ["synth_val_c"]


def test_trial_log_records_every_run(planted_bench, capsys):
    b = planted_bench
    mod = _load_script("bench_detector.py")
    params = {"log": {"min_snr": 2.5}}
    pfile = b["root"] / "params.json"
    pfile.write_text(json.dumps(params), encoding="utf-8")
    assert mod.main(["--backend", "log", "--params", str(pfile), "--split", "train", "--tags", "synth_HS_b",
                     "--out", str(b["root"] / "o1"), *_bench_args(b)]) == 0
    assert mod.main(["--backend", "log", "--params", json.dumps(params), "--split", "val",
                     "--out", str(b["root"] / "o2"), *_bench_args(b)]) == 0
    lines = bd.read_jsonl(b["root"] / "trials.jsonl")
    assert len(lines) == 2
    first = lines[0]
    for k in ("utc", "git_sha", "backend", "params", "params_sha256", "split", "tags", "pooled", "per_gonad",
              "skipped", "final", "bootstrap", "versions"):
        assert k in first, k
    assert first["backend"] == "log" and first["params"] == params and first["split"] == "train"
    assert first["params_sha256"] == hashlib.sha256(json.dumps(params, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert first["tags"] == ["synth_HS_b"] and first["final"] is False
    assert first["utc"].endswith("+00:00") and "f1" in first["pooled"]["spots_all"]
    assert lines[1]["split"] == "val" and lines[1]["tags"] == ["synth_val_c"]
    assert not (b["root"] / "final_runs.jsonl").exists()
    # --no-trial-log leaves the log alone
    assert mod.main(["--backend", "log", "--split", "val", "--no-trial-log", "--out", str(b["root"] / "o3"),
                     *_bench_args(b)]) == 0
    assert len(bd.read_jsonl(b["root"] / "trials.jsonl")) == 2


def test_guard_refuses_test_gonads(planted_bench, capsys):
    b = planted_bench
    mod = _load_script("bench_detector.py")
    out = b["root"] / "o"
    assert mod.main(["--backend", "log", "--split", "test", "--out", str(out), *_bench_args(b)]) == 2
    assert "REFUSED" in capsys.readouterr().out
    # a test tag named under another split is refused as being outside it; under test, without --final
    assert mod.main(["--backend", "log", "--split", "train", "--tags", "synth_test_d", "--out", str(out),
                     *_bench_args(b)]) == 2
    assert mod.main(["--backend", "log", "--split", "test", "--tags", "synth_test_d", "--out", str(out),
                     *_bench_args(b)]) == 2
    assert mod.main(["--backend", "log", "--split", "test", "--final", "--out", str(out), *_bench_args(b)]) == 2
    assert "--reason" in capsys.readouterr().out
    # --final can never skip the final log
    assert mod.main(["--backend", "log", "--split", "test", "--final", "--reason", "x", "--allow-dirty",
                     "--no-trial-log", "--out", str(out), *_bench_args(b)]) == 2
    assert "--no-trial-log is not allowed with --final" in capsys.readouterr().out
    assert not out.exists() and not (b["root"] / "trials.jsonl").exists()
    # the library refuses as well, before it loads anything
    with pytest.raises(bd.TestSetGuardError):
        bd.score("synth_test_d", pd.DataFrame(), pd.DataFrame(), manifest=b["manifest"], cache_root=b["root"] / "cache")
    with pytest.raises(bd.TestSetGuardError):
        bd.run_backend("n2rad51_HS_male_001", "log", {}, cache_root=b["root"] / "nowhere")
    with pytest.raises(bd.TestSetGuardError):
        bd.guard_test(["n2rad51_noHS_herm_000", "n2rad51_noHS_male_002"])
    bd.guard_test(["n2rad51_noHS_herm_000", "n2rad51_HS_herm_001"])          # train and val pass
    assert bd.is_test_tag("synth_x", "test") and not bd.is_test_tag("synth_x", "val")
    with pytest.raises(ValueError, match="pre-registered"):
        bd.gonad_role("n2rad51_HS_male_001", "train")                          # a manifest cannot move a gonad
    with pytest.raises(ValueError, match="unusable"):
        bd.score("synth_rg_e", pd.DataFrame(), pd.DataFrame(), manifest=b["manifest"], cache_root=b["root"] / "cache")


def test_final_runs_go_to_their_own_log_and_are_limited(planted_bench, capsys):
    """--final on a SYNTHETIC test gonad (temporary logs): the record lands in final_runs.jsonl, not in
    trials.jsonl, and a third test scoring of the same backend is refused."""
    b = planted_bench
    mod = _load_script("bench_detector.py")
    args = ["--backend", "log", "--split", "test", "--final", "--reason", "synthetic check of the final path",
            "--allow-dirty", *_bench_args(b)]
    assert mod.main([*args, "--out", str(b["root"] / "f1")]) == 0
    finals = bd.read_jsonl(b["root"] / "final_runs.jsonl")
    assert len(finals) == 1 and finals[0]["final"] is True and finals[0]["reason"].startswith("synthetic")
    assert finals[0]["tags"] == ["synth_test_d"] and not (b["root"] / "trials.jsonl").exists()
    assert mod.main([*args, "--out", str(b["root"] / "f2")]) == 0
    capsys.readouterr()
    assert mod.main([*args, "--out", str(b["root"] / "f3")]) == 2
    assert "allows 2" in capsys.readouterr().out
    assert len(bd.read_jsonl(b["root"] / "final_runs.jsonl")) == 2


def test_runner_process_pool(planted_bench, capsys):
    b = planted_bench
    mod = _load_script("bench_detector.py")
    out = b["root"] / "pool"
    assert mod.main(["--backend", "log", "--split", "train", "--jobs", "2", "--no-trial-log", "--out", str(out),
                     *_bench_args(b)]) == 0
    seq = b["root"] / "seq"
    assert mod.main(["--backend", "log", "--split", "train", "--no-trial-log", "--out", str(seq), *_bench_args(b)]) == 0
    a = pd.read_csv(out / "per_gonad.csv").drop(columns=["detect_s", "score_s", "elapsed_s"], errors="ignore")
    s = pd.read_csv(seq / "per_gonad.csv").drop(columns=["detect_s", "score_s", "elapsed_s"], errors="ignore")
    pd.testing.assert_frame_equal(a, s)


def test_cache_skips_entries_without_outputs(tmp_path):
    e = Entry(tag="synth_empty", results_dir=tmp_path / "nothing", split="train", nd2=tmp_path / "x.nd2",
              stem="x", extra={"role": "train"})
    st = bd.cache_gonad(e, tmp_path / "cache")
    assert st["status"] == "skipped" and "no pipeline outputs" in st["reason"]
    with pytest.raises(FileNotFoundError, match="no complete cache"):
        bd.load_cache("synth_empty", tmp_path / "cache")
    assert bd.cache_size_bytes(tmp_path / "cache") == {"per_tag": {}, "total": 0}
