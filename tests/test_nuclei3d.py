"""The EXPERIMENTAL nuclei3d nucleus model (germquant.segment.nuclei3d) and the protocol scope of the nucleus
benchmark, on the CPU with tiny networks and synthetic volumes:

- the network's shapes, heads and receptive radius (a perturbation never reaches farther);
- the loss: weights, the zero-weight voxels, the normalisation;
- the targets on synthetic masks: cell distance, neighbour closeness, the class map of the partial
  annotation (traced, shell, ignore inside untraced Cellpose labels plus the margin, background);
- the instance step: h-maxima seeds and the seeded watershed split two touching spheres;
- tiled inference equals untiled inference (xy and z);
- the backend contract (segment_nuclei, the provenance, the pipeline's sub-hash and manifest, the default
  path untouched);
- ``germquant train-nuclei3d`` for one epoch on a tiny synthetic bench cache, never opening the test gonad;
- the protocol mode of germquant.benchmark.nuclei3d and scripts/benchmark_nuclei3d.py (the scope rule,
  fragments, the 20 percent splits and merges, the test guard, --pred-dir).
Nothing here needs microscopy data or a GPU.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from germquant.benchmark import nuclei3d as n3  # noqa: E402
from germquant.segment.nuclei3d import data as d3  # noqa: E402
from germquant.segment.nuclei3d import infer as inf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SP = (0.2, 0.108, 0.108)


def _model(levels=3, base=4, seed=0):
    from germquant.segment.nuclei3d.model import build_model

    torch.manual_seed(seed)
    m = build_model({"base_channels": base, "levels": levels})
    g = torch.Generator().manual_seed(seed + 1)
    for mod in m.modules():
        if isinstance(mod, torch.nn.BatchNorm3d):
            mod.running_mean.copy_(torch.randn(mod.num_features, generator=g) * 0.1)
            mod.running_var.copy_(torch.rand(mod.num_features, generator=g) + 0.5)
    return m.eval()


def _spheres(shape=(24, 64, 64), centres=((12, 32, 20), (12, 32, 44)), r_um=1.3, sp=SP):
    """Label image of spheres (radius in um) at voxel centres, the earlier label winning on overlap."""
    zz, yy, xx = np.indices(shape).astype(float)
    lab = np.zeros(shape, np.int32)
    for k, c in enumerate(centres, start=1):
        d2 = ((zz - c[0]) * sp[0]) ** 2 + ((yy - c[1]) * sp[1]) ** 2 + ((xx - c[2]) * sp[2]) ** 2
        lab[(d2 <= r_um ** 2) & (lab == 0)] = k
    return lab


def _save_random_model(path, levels=3, base=4, post=None):
    m = _model(levels=levels, base=base)
    meta = {"input_roles": ["dna"], "norm": {"low_pct": 1.0, "high_pct": 99.9, "stride_yx": 4},
            "post": post or {"fg_threshold": 0.5, "h_um": 0.4, "neighbour_weight": 1.0},
            "spacing_zyx_um": list(SP)}
    return inf.save_model_file(path, m, meta)


# ------------------------------------------------------------------------------------ model
def test_model_shapes_heads_and_radius():
    from germquant.segment.nuclei3d.model import HEAD_BIAS, HEADS, receptive_radius

    m = _model(levels=4, base=4)
    out = m(torch.zeros(2, 1, 8, 32, 24))
    assert tuple(out.shape) == (2, len(HEADS), 8, 32, 24) and HEADS[0] == "foreground"
    assert m.head.bias.detach().tolist() == list(HEAD_BIAS)
    assert m.arch() == {"in_channels": 1, "base_channels": 4, "levels": 4}
    # a perturbed input voxel never changes an output farther away than the receptive radius
    m3 = _model(levels=3, base=4)
    shape = (32, 64, 64)
    x = torch.zeros((1, 1) + shape)
    c = tuple(s // 2 for s in shape)
    x2 = x.clone()
    x2[(0, 0) + c] = 10.0
    with torch.no_grad():
        d = (m3(x2) - m3(x)).abs()[0].amax(dim=0).numpy()
    hit = np.argwhere(d > 0)
    assert len(hit) > 0
    assert np.all(np.abs(hit - np.array(c)).max(axis=0) <= np.array(receptive_radius(3)))


def test_loss_weights_and_normalisation():
    from germquant.segment.nuclei3d.train import nuclei3d_loss

    g = torch.Generator().manual_seed(0)
    out = torch.randn(2, 3, 4, 6, 6, generator=g)
    tgt = torch.zeros(2, 3, 4, 6, 6)
    tgt[:, 0, 1:3, 2:4, 2:4] = 1.0
    tgt[:, 1] = tgt[:, 0] * 0.7
    tgt[:, 2] = tgt[:, 0] * 0.3
    w = torch.ones(2, 1, 4, 6, 6)
    w[:, :, 0] = 0.0                                    # one plane ignored
    loss, terms = nuclei3d_loss(out, tgt, w, lambda_dist=2.0, lambda_nb=0.5)
    bce = torch.nn.functional.binary_cross_entropy_with_logits(out[:, 0], tgt[:, 0], reduction="none")
    manual = ((bce + 2.0 * (out[:, 1] - tgt[:, 1]) ** 2 + 0.5 * (out[:, 2] - tgt[:, 2]) ** 2) * w[:, 0]).sum() / w.sum()
    assert torch.allclose(loss, manual, atol=1e-6)
    assert abs(float(terms["bce"]) - float((bce * w[:, 0]).sum() / w.sum())) < 1e-6
    # whatever the network says where the weight is 0 changes nothing
    out2 = out.clone()
    out2[:, :, 0] = 1000.0
    assert torch.allclose(nuclei3d_loss(out2, tgt, w, lambda_dist=2.0, lambda_nb=0.5)[0], loss, atol=1e-6)
    # weights of shape (B, Z, Y, X) work too; an all-zero weight gives 0, not a division by zero
    assert torch.allclose(nuclei3d_loss(out, tgt, w[:, 0], lambda_dist=2.0, lambda_nb=0.5)[0], loss)
    assert float(nuclei3d_loss(out, tgt, torch.zeros_like(w))[0]) == 0.0


# ---------------------------------------------------------------------------------- targets
def test_targets_on_synthetic_masks():
    from scipy import ndimage as ndi

    truth = _spheres()
    # Cellpose: the two traced nuclei, each a little larger (radius 1.5 um), and an untraced third nucleus
    cp = _spheres(centres=((12, 32, 20), (12, 32, 44), (12, 8, 32)), r_um=1.5)
    t = d3.build_targets(truth, cp, SP, {})
    cls, dist, nb = t["cls"], t["dist"], t["nb"]
    info = t["info"]
    assert info["n_traced"] == 2 and info["n_cellpose"] == 3 and info["n_cellpose_traced"] == 2
    # classes
    assert np.all(cls[truth > 0] == d3.CLASS_TRACED)
    assert np.all(cls[cp == 3] == d3.CLASS_IGNORE)
    shell = (cp > 0) & (cp < 3) & (truth == 0)
    margin = d3.dilate_ellipsoid(cp == 3, 0.5, SP)
    assert np.all(cls[shell & ~margin] == d3.CLASS_SHELL) and (shell & ~margin).any()
    assert np.all(cls[margin & (truth == 0)] == d3.CLASS_IGNORE)
    assert np.all(cls[(cp == 0) & ~margin] == d3.CLASS_BACKGROUND)
    w = d3.class_weights({})[cls]
    assert w[cls == d3.CLASS_IGNORE].max() == 0 and w[cls == d3.CLASS_SHELL].min() == 1.0
    assert d3.class_weights({"w_traced_shell": 0.0})[d3.CLASS_SHELL] == 0.0
    # cell distance: the EDT inside each sphere (with the other sphere and the stack edge as outside)
    for k in (1, 2):
        ref = ndi.distance_transform_edt(np.pad(truth == k, 1), sampling=SP)[1:-1, 1:-1, 1:-1]
        np.testing.assert_allclose(dist[truth == k], ref[truth == k], rtol=0, atol=1e-6)
    assert np.all(dist[truth == 0] == 0) and 1.1 < dist.max() < 1.35
    # neighbour closeness: high at the contact between the two spheres, zero far from any neighbour
    contact = (truth == 1) & ndi.binary_dilation(truth == 2)
    assert contact.any() and nb[contact].min() > 0.5
    # (the neighbours of nucleus 1 are the other traced mask and every Cellpose label that is not its own,
    # so the larger Cellpose label of nucleus 2 counts with its full extent)
    others = (truth == 2) | (cp == 2) | (cp == 3)
    far = (truth == 1) & (ndi.distance_transform_edt(~others, sampling=SP) > 1.05)
    assert far.any() and np.all(nb[far] == 0)
    near = (truth == 1) & (ndi.distance_transform_edt(~others, sampling=SP) < 0.95)
    assert near.any() and np.all(nb[near] > 0)
    assert np.all(nb[truth == 0] == 0) and nb.max() <= 1.0
    # a traced nucleus's own (larger) Cellpose label is not its neighbour: alone, it has no closeness at all
    one = _spheres(centres=((12, 32, 20),))
    t1 = d3.build_targets(one, _spheres(centres=((12, 32, 20),), r_um=1.6), SP, {})
    assert t1["nb"].max() == 0 and t1["info"]["n_cellpose_traced"] == 1


def test_targets_cache_round_trip(tmp_path):
    import tifffile

    truth = _spheres()
    tp = tmp_path / "t__imaris_labels.tif"
    tifffile.imwrite(tp, truth)
    meta = {"arrays": {"labels": {"data_sha256": "abc"}}}
    a = d3.load_targets(tmp_path / "prep", "g", meta, truth, tp, SP, {})
    b = d3.load_targets(tmp_path / "prep", "g", meta, truth, tp, SP, {})       # from the files
    assert isinstance(b["cls"], np.memmap)
    np.testing.assert_array_equal(np.asarray(a["cls"]), np.asarray(b["cls"]))
    np.testing.assert_array_equal(np.asarray(a["dist"]), np.asarray(b["dist"]))
    assert b["traced_centers"].shape == (2, 3)
    import gc

    del a, b                     # Windows cannot replace a file that is still memory-mapped
    gc.collect()
    c = d3.load_targets(tmp_path / "prep", "g", meta, truth, tp, SP, {"margin_um": 0.3})   # stale: rebuilt
    assert c["info"]["settings"]["margin_um"] == 0.3


def test_patches_are_reproducible_and_weighted():
    truth = _spheres(shape=(24, 96, 96), centres=((12, 48, 36), (12, 48, 60)))
    t = d3.build_targets(truth, truth, SP, {})
    dna = (100 + 600 * (truth > 0)).astype(np.uint16)
    g = d3.TrainGonad(tag="g", spacing=SP, dna=dna, lo_hi=inf.robust_range(dna), cls=t["cls"],
                      dist=t["dist"].astype(np.float16), nb=t["nb"].astype(np.float16),
                      traced_centers=d3.label_centroids(truth)[2], germline_centers=np.zeros((0, 3)))
    a = d3.sample([g], (8, 32, 32), {}, seed=1, epoch=2, index=3)
    b = d3.sample([g], (8, 32, 32), {}, seed=1, epoch=2, index=3)
    for u, v in zip(a, b):
        np.testing.assert_array_equal(u, v)
    x, tgt, w = d3.make_patch(g, (8, 30, 20), (8, 32, 32), {}, None, augment=False)
    assert x.shape == (1, 8, 32, 32) and tgt.shape == (3, 8, 32, 32) and w.shape == (1, 8, 32, 32)
    np.testing.assert_array_equal(tgt[0], (truth[8:16, 30:62, 20:52] > 0).astype(np.float32))
    # hanging over the edge: zero input and zero weight outside the stack
    x2, tgt2, w2 = d3.make_patch(g, (-4, -8, 0), (8, 32, 32), {}, None, augment=False)
    assert np.all(w2[:, :4] == 0) and np.all(x2[:, :4] == 0) and np.all(w2[:, :, :8, :] == 0)
    assert w2[:, 4:, 8:, :].any()
    batches = list(d3.Prefetch(d3.batches([g], (8, 32, 32), {}, 0, 0, 3, 2)))
    assert len(batches) == 3 and batches[0][1].shape == (2, 3, 8, 32, 32)


# ---------------------------------------------------------------------------- the instances
def test_seeded_watershed_splits_two_touching_spheres():
    truth = _spheres()
    t = d3.build_targets(truth, truth, SP, {})
    maps = {"fg": (truth > 0).astype(np.float32), "dist": t["dist"], "nb": t["nb"]}
    lab = inf.instances(maps, h_um=0.4)
    assert lab.dtype == np.int32 and lab.max() == 2
    for k in (1, 2):                                # each instance is one sphere, nearly voxel for voxel
        hit = lab[truth == k]
        best = np.bincount(hit[hit > 0]).argmax()
        assert (hit == best).mean() > 0.97
    # the distance map alone (no neighbour damping) still splits them; one sphere stays one
    assert inf.instances(maps, h_um=0.4, neighbour_weight=0.0).max() == 2
    one = _spheres(centres=((12, 32, 32),))
    t1 = d3.build_targets(one, one, SP, {})
    assert inf.instances({"fg": (one > 0).astype(np.float32), "dist": t1["dist"], "nb": t1["nb"]}).max() == 1
    # an h above every dome leaves no seed: nothing is labelled
    assert inf.instances(maps, h_um=5.0).max() == 0
    # the foreground threshold applies; an empty foreground gives an empty image
    assert inf.instances({**maps, "fg": maps["fg"] * 0.4}, fg_threshold=0.5).max() == 0
    # deterministic
    np.testing.assert_array_equal(inf.instances(maps), inf.instances(maps))


def test_tiled_inference_equals_untiled_in_xy_and_z():
    rng = np.random.default_rng(3)
    m = _model(levels=3, base=4)
    dna = rng.integers(90, 900, size=(20, 100, 90)).astype(np.uint16)
    lohi = inf.robust_range(dna)
    whole = inf.predict_maps(m, dna, lohi, tile_zyx=None)
    tiled = inf.predict_maps(m, dna, lohi, tile_zyx=(96, 64, 64))
    for k in ("fg", "dist", "nb"):
        assert whole[k].shape == dna.shape and whole[k].dtype == np.float32
        np.testing.assert_allclose(tiled[k], whole[k], rtol=0, atol=2e-6)
    assert whole["dist"].min() >= 0 and 0 <= whole["nb"].min() and whole["nb"].max() <= 1
    d2 = rng.integers(90, 900, size=(40, 24, 24)).astype(np.uint16)
    r2 = inf.robust_range(d2)
    a = inf.predict_maps(m, d2, r2, tile_zyx=(32, None, None))
    b = inf.predict_maps(m, d2, r2, tile_zyx=None)
    for k in a:
        np.testing.assert_allclose(a[k], b[k], rtol=0, atol=2e-6)
    # the seed-map path equals the maps path
    fg, s = inf.predict_seed_map(m, dna, lohi, tile_zyx=(96, 64, 64))
    fg2, s2 = inf.seed_map_block(whole["fg"], whole["dist"], whole["nb"], 0.5, 1.0)
    np.testing.assert_array_equal(fg, fg2)
    np.testing.assert_allclose(s, s2, rtol=0, atol=2e-6)


# ------------------------------------------------------------------------------- backend
def test_backend_contract_with_a_random_model(tmp_path, caplog):
    from germquant import provenance
    from germquant.segment import backend_provenance, list_backends, segment_nuclei

    assert "nuclei3d" in list_backends()
    path = tmp_path / "model.pt"
    sha = _save_random_model(path, post={"fg_threshold": 0.0, "h_um": 0.05, "neighbour_weight": 1.0})
    assert inf.read_model_file(path)["input_roles"] == ["dna"]
    truth = _spheres(shape=(16, 64, 64), centres=((8, 32, 20), (8, 32, 44)))
    dna = (100 + 700 * (truth > 0)).astype(np.uint16)
    section = {"nuclei3d": {"model": str(path), "device": "cpu", "tile_zyx": [96, 64, 64]}}
    labels, used = segment_nuclei(dna, SP, method="cellpose", backend="nuclei3d", min_volume_um3=0.01,
                                  cfg_section=section)
    assert used == "nuclei3d" and labels.shape == dna.shape and labels.dtype == np.int32
    again, _ = segment_nuclei(dna, SP, method="cellpose", backend="nuclei3d", min_volume_um3=0.01, cfg_section=section)
    np.testing.assert_array_equal(labels, again)
    # the config keys override the model file's instance parameters: a threshold above every probability
    none, _ = segment_nuclei(dna, SP, method="cellpose", backend="nuclei3d", min_volume_um3=0.01,
                             cfg_section={"nuclei3d": {**section["nuclei3d"], "fg_threshold": 1.0}})
    assert none.max() == 0
    # no model: an error under method cellpose, the classical fallback under auto
    with pytest.raises(ValueError, match="nuclei3d.model"):
        segment_nuclei(dna, SP, method="cellpose", backend="nuclei3d", cfg_section={})
    _, fallback = segment_nuclei(dna.astype(np.float32), SP, method="auto", backend="nuclei3d", cfg_section={},
                                 min_volume_um3=0.01)
    assert fallback == "classical"
    with pytest.raises(FileNotFoundError):
        segment_nuclei(dna, SP, method="cellpose", backend="nuclei3d",
                       cfg_section={"nuclei3d": {"model": str(tmp_path / "missing.pt"), "device": "cpu"}})
    # provenance: the model file's sha256, never fatal; nothing for the built-in default backends
    rec = backend_provenance("nuclei3d", section)
    assert rec == {"backend": "nuclei3d", "model": str(path), "model_sha256": sha} and sha == provenance.file_sha256(path)
    assert backend_provenance("nuclei3d", {})["model_sha256"] == "absent"
    assert backend_provenance("cellpose", section) == {} and backend_provenance("classical", {}) == {}
    assert backend_provenance("no_such_backend", {}) == {}
    caplog.set_level("WARNING")
    segment_nuclei(dna, SP, method="cellpose", backend="nuclei3d", min_volume_um3=0.01,
                   cfg_section={"nuclei3d": {**section["nuclei3d"], "bogus": 1}})
    assert any("unknown keys ['bogus']" in r.getMessage() for r in caplog.records)
    # a file that is not a nuclei3d model is refused
    from germquant.spots.net3d.infer import save_model_file as save_net3d
    from germquant.spots.net3d.model import build_model as build_net3d

    save_net3d(tmp_path / "net3d.pt", build_net3d({"base_channels": 4, "levels": 3}), {})
    with pytest.raises(ValueError, match="not a nuclei3d model"):
        inf.read_model_file(tmp_path / "net3d.pt")


def test_default_path_hashes_and_manifest_are_unchanged(tmp_path):
    from germquant import pipeline
    from germquant.config import load_config
    from germquant.stages import stage_hashes

    cfg = load_config(ROOT / "config" / "config.yaml")
    roles = {"dna": 0, "central_element": 1, "foci": 2, "granule": None}
    assert pipeline.segment_backend_provenance(cfg) == {}
    tools = {"cellpose": "x"}
    assert stage_hashes(cfg, roles, tools, "m") == stage_hashes(cfg, roles, tools, "m", segment_backend={})
    path = tmp_path / "m.pt"
    sha = _save_random_model(path)
    cfg.set("segmentation.nuclei.backend", "nuclei3d")
    cfg.set("segmentation.nuclei.nuclei3d", {"model": str(path)})
    prov = pipeline.segment_backend_provenance(cfg)
    assert prov == {"backend": "nuclei3d", "model": str(path), "model_sha256": sha}
    h1 = stage_hashes(cfg, roles, tools, "m", segment_backend=prov)
    h0 = stage_hashes(cfg, roles, tools, "m")
    assert h1["segment"] != h0["segment"]
    assert {k: v for k, v in h1.items() if k != "segment"} == {k: v for k, v in h0.items() if k != "segment"}
    # method classical never calls the backend: no provenance
    cfg.set("segmentation.nuclei.method", "classical")
    assert pipeline.segment_backend_provenance(cfg) == {}


def test_pipeline_runs_nuclei3d(tmp_path, monkeypatch):
    from test_end_to_end import CHANNELS, SPACING, _synthetic_stack

    from germquant import pipeline
    from germquant.config import load_config

    stack = _synthetic_stack()
    meta = {"sizes": {"C": 3, "Z": 20, "Y": 80, "X": 320}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)
    path = tmp_path / "m.pt"
    sha = _save_random_model(path, post={"fg_threshold": 0.0, "h_um": 0.05, "neighbour_weight": 1.0})
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg.set("segmentation.nuclei.method", "cellpose")
    cfg.set("segmentation.nuclei.cellpose_model", "cpsam")
    cfg.set("segmentation.nuclei.backend", "nuclei3d")
    cfg.set("segmentation.nuclei.nuclei3d", {"model": str(path), "device": "cpu", "tile_zyx": [96, 128, 128]})
    cfg.set("spots.enabled", False)
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["segment"]["status"] == "ran", res["stages"]["segment"]
    isum = res["tables"]["image_summary"]
    assert isum["segmentation_method"].iloc[0] == "nuclei3d"
    man = json.loads((tmp_path / "r" / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["segment_backend"] == {"backend": "nuclei3d", "model": str(path), "model_sha256": sha}


# ------------------------------------------------------------------- the CLI and the harness
@pytest.fixture
def tiny_bench(tmp_path):
    from test_bench_detector import _write_gonad

    root = tmp_path
    gonads = [_write_gonad(root, "synth_noHS_a", "train", "noHS", (2, 4, 6, 3), 11),
              _write_gonad(root, "synth_HS_b", "train", "HS", (5, 1, 3, 7), 12),
              _write_gonad(root, "synth_val_c", "val", "noHS", (3, 3, 5, 2), 13),
              _write_gonad(root, "synth_test_d", "test", "HS", (4, 2, 2, 6), 14)]
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({"gonads": [g["entry"] for g in gonads]}), encoding="utf-8")
    split = root / "split.json"
    split.write_text(json.dumps({"roles": {"train": ["synth_noHS_a", "synth_HS_b"], "val": ["synth_val_c"],
                                           "test": ["synth_test_d"]}}), encoding="utf-8")
    return {"root": root, "manifest": manifest, "split": split}


def _cli_args(b, out, *extra):
    return ["train-nuclei3d", "--cache", str(b["root"] / "cache"), "--split-file", str(b["split"]),
            "--out", str(out), "--epochs", "1", "--patch", "16,32,32", "--seed", "3",
            "--manifest", str(b["manifest"]), "--batch-size", "2", "--steps-per-epoch", "2", "--device", "cpu",
            "--base-channels", "4", "--levels", "3", "--tile", "96,128,128", *extra]


def test_train_nuclei3d_cli_one_epoch_on_a_tiny_cache(tiny_bench, monkeypatch, capsys):
    from germquant import cli
    from germquant.benchmark import detector as bd
    from germquant.segment.nuclei3d.train import resolve_gonads

    b = tiny_bench
    opened = []
    real_load = bd.load_cache

    def spy_load(tag, *a, **k):
        opened.append(tag)
        return real_load(tag, *a, **k)

    monkeypatch.setattr(bd, "load_cache", spy_load)
    out = b["root"] / "model"
    # after two steps the network is near its initial biases (foreground probability about 0.12), so the
    # grid includes a low foreground threshold that gives the scorer something to match
    rc = cli.main(_cli_args(b, out, "--post-grid", "fg_threshold=0.05,0.5;h_um=0.05,0.3"))
    text = capsys.readouterr().out
    assert rc == 0, text
    assert "synth_test_d" not in opened and set(opened) == {"synth_noHS_a", "synth_HS_b", "synth_val_c"}
    for f in ("model.pt", "training_log.json", "model_card.json", "checkpoints/epoch_001.pt"):
        assert (out / f).is_file(), f
    logd = json.loads((out / "training_log.json").read_text(encoding="utf-8"))
    card = json.loads((out / "model_card.json").read_text(encoding="utf-8"))
    assert logd["status"] == "finished" and logd["best_epoch"] == 1 and len(logd["epochs"]) == 1
    assert logd["model_sha256"] == card["model_sha256"] == inf.file_sha256(out / "model.pt")
    assert logd["checkpoints"][0]["sha256"] == inf.file_sha256(out / "checkpoints" / "epoch_001.pt")
    ev = logd["epochs"][0]["val"]
    assert set(ev["per_gonad"]) == {"synth_val_c"} and "f1_iou50" in ev["pooled"]
    assert ev["pooled"]["n_truth"] == 4                         # the four planted nuclei of the val gonad
    assert logd["post_grid"]["n_configurations"] == 4
    assert logd["selection"]["n_configurations_scored_on_val"] == 5
    assert [g["tag"] for g in card["training_gonads"]] == ["synth_noHS_a", "synth_HS_b"]
    assert card["test_gonads_never_opened"] == ["synth_test_d"] and card["git_sha"]
    doc = inf.read_model_file(out / "model.pt")
    assert doc["post"] == logd["selection"]["post"] == logd["post_grid"]["chosen"]
    assert doc["post"]["h_um"] in (0.05, 0.3) and doc["post"]["fg_threshold"] in (0.05, 0.5)
    assert doc["input_roles"] == ["dna"] and doc["training"]["train_gonads"] == ["synth_noHS_a", "synth_HS_b"]
    # the same seed gives the same weights
    out2 = b["root"] / "model2"
    assert cli.main(_cli_args(b, out2, "--prep-dir", str(out / "prep"))) == 0
    w1 = inf.read_model_file(out / "checkpoints" / "epoch_001.pt")["state_dict"]
    w2 = inf.read_model_file(out2 / "checkpoints" / "epoch_001.pt")["state_dict"]
    assert all(torch.equal(w1[k], w2[k]) for k in w1)
    # a split file that moves a registered gonad, or disagrees with the manifest, is refused
    bad = b["root"] / "bad.json"
    bad.write_text(json.dumps({"roles": {"train": ["synth_noHS_a"], "val": ["synth_HS_b"], "test": []}}),
                   encoding="utf-8")
    with pytest.raises(ValueError, match="manifest says"):
        resolve_gonads(bad, b["manifest"])
    moved = b["root"] / "moved.json"
    moved.write_text(json.dumps({"roles": {"train": ["n2rad51_HS_male_001"]}}), encoding="utf-8")
    with pytest.raises(ValueError, match="pre-registered"):
        resolve_gonads(moved, b["manifest"])

    # the benchmark script scores the model on val through the pipeline's path, then again from --pred-dir
    spec = importlib.util.spec_from_file_location("benchmark_nuclei3d", ROOT / "scripts" / "benchmark_nuclei3d.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    bench = b["root"] / "bench"
    rc = mod.main(["--manifest", str(b["manifest"]), "--out", str(bench), "--coverage", "protocol",
                   "--model", str(out / "model.pt"), "--cache-root", str(b["root"] / "cache"), "--device", "cpu",
                   "--tile", "96,128,128", "--n-boot", "20"])
    assert rc == 0, capsys.readouterr().out
    pg = pd.read_csv(bench / "per_gonad.csv")
    pooled = json.loads((bench / "pooled.json").read_text(encoding="utf-8"))
    # train gonads are scored if asked; the test gonad never is
    assert "synth_test_d" not in pg["tag"].tolist()
    assert [s["tag"] for s in pooled["skipped"]] == ["synth_test_d"]
    assert pg["coverage_rule"].str.startswith("protocol").all() and pg["pred_source"].eq("nuclei3d model").all()
    assert (bench / "labels" / "synth_val_c__labels.tif").is_file()
    rc = mod.main(["--manifest", str(b["manifest"]), "--out", str(b["root"] / "bench2"), "--coverage", "protocol",
                   "--pred-dir", str(bench / "labels"), "--tags", "synth_val_c", "--n-boot", "20"])
    assert rc == 0
    pg2 = pd.read_csv(b["root"] / "bench2" / "per_gonad.csv")
    row1 = pg[pg["tag"] == "synth_val_c"].iloc[0]
    np.testing.assert_equal(pg2["f1_iou50"].iloc[0], row1["f1_iou50"])
    assert pg2["tp_iou50"].iloc[0] == row1["tp_iou50"] and pg2["n_pred"].iloc[0] == row1["n_pred"]
    # the validation inside training scored the same labels the same way
    chosen = next(r for r in logd["post_grid"]["results"] if r["post"] == logd["post_grid"]["chosen"])
    for k in ("tp_iou50", "fp_iou50", "fn_iou50", "n_pred", "n_pred_counted"):
        assert chosen["per_gonad"]["synth_val_c"][k] == int(row1[k]), k


def _oracle_maps(dna):
    """Ideal maps of a planted bench gonad (DAPI 1000 inside the nuclei, 50 outside): the foreground, and
    the cell distance and neighbour closeness of its connected components."""
    from scipy import ndimage as ndi

    lab, _ = ndi.label(np.asarray(dna) > 500)
    t = d3.build_targets(lab, lab, SP, {})
    return {"fg": (lab > 0).astype(np.float32), "dist": t["dist"], "nb": t["nb"]}


def test_validation_scorer_equals_the_benchmark_script(tiny_bench, monkeypatch):
    """With an oracle in place of the network, the training-time validation (evaluate, post_grid) and
    scripts/benchmark_nuclei3d.py --model score the same labels identically, and perfectly."""
    from germquant.segment.nuclei3d import train as tr

    def oracle_seed_map(model, dna, lo_hi, *, device="cpu", tile_zyx=None, fg_threshold=0.5, neighbour_weight=1.0):
        m = _oracle_maps(dna)
        return inf.seed_map_block(m["fg"], m["dist"], m["nb"], fg_threshold, neighbour_weight)

    def oracle_maps(model, dna, lo_hi, *, device="cpu", tile_zyx=None):
        return _oracle_maps(dna)

    monkeypatch.setattr(inf, "predict_seed_map", oracle_seed_map)
    monkeypatch.setattr(tr, "predict_seed_map", oracle_seed_map)
    monkeypatch.setattr(tr, "predict_maps", oracle_maps)
    b = tiny_bench
    gon = tr.resolve_gonads(b["split"], b["manifest"])
    vals = [tr.load_val_gonad(t, gon["entries"][t], b["root"] / "cache", norm=dict(inf.NORM_DEFAULTS))
            for t in gon["val"]]
    ev = tr.evaluate(None, vals)
    pv = ev["pooled"]
    assert pv["tp_iou50"] == 4 and pv["fp_iou50"] == 0 and pv["fn_iou50"] == 0 and pv["f1_iou50"] == 1.0
    assert pv["ap_iou50_95"] == 1.0 and pv["n_splits"] == 0 and pv["n_merges"] == 0
    grid = tr.post_grid(None, vals, {"h_um": [0.4, 5.0]})
    assert grid[0]["pooled"]["tp_iou50"] == 4 and grid[1]["pooled"]["tp_iou50"] == 0
    assert grid[0]["per_gonad"] == ev["per_gonad"]
    # the script, through segment_nuclei and the nuclei3d backend, on the same oracle
    path = b["root"] / "m.pt"
    _save_random_model(path)
    spec = importlib.util.spec_from_file_location("benchmark_nuclei3d", ROOT / "scripts" / "benchmark_nuclei3d.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    out = b["root"] / "bench_oracle"
    rc = mod.main(["--manifest", str(b["manifest"]), "--out", str(out), "--coverage", "protocol", "--split", "val",
                   "--model", str(path), "--cache-root", str(b["root"] / "cache"), "--device", "cpu",
                   "--n-boot", "20"])
    assert rc == 0
    row = pd.read_csv(out / "per_gonad.csv").iloc[0]
    for k in ("tp_iou50", "fp_iou50", "fn_iou50", "n_pred", "n_pred_counted", "n_splits", "n_merges"):
        assert int(row[k]) == ev["per_gonad"]["synth_val_c"][k], k
    assert abs(row["ap_iou50_95"] - ev["per_gonad"]["synth_val_c"]["ap_iou50_95"]) < 1e-12


# ----------------------------------------------------------------------- protocol scoring
def _cube(vol, label, z, y, x):
    vol[z[0]:z[1], y[0]:y[1], x[0]:x[1]] = label


def test_protocol_scope_rule_fragments_and_twenty_percent_rule():
    vox = (0.2, 0.1, 0.1)
    truth = np.zeros((10, 60, 60), np.int32)
    _cube(truth, 1, (2, 6), (5, 15), (5, 15))           # 400 voxels
    _cube(truth, 2, (2, 6), (5, 15), (15, 25))          # 400 voxels, touching 1
    _cube(truth, 3, (0, 1), (50, 52), (50, 52))         # a 4-voxel fragment
    pred = np.zeros_like(truth)
    _cube(pred, 10, (2, 6), (5, 15), (5, 12))           # 280 voxels inside truth 1 (70 percent of it)
    _cube(pred, 11, (2, 6), (5, 15), (12, 15))          # 120 voxels inside truth 1 (30 percent)
    _cube(pred, 12, (2, 6), (5, 15), (15, 25))          # truth 2 exactly
    _cube(pred, 13, (2, 6), (20, 26), (5, 15))          # touching C only partly: rows 20 to 25
    _cube(pred, 14, (2, 6), (40, 48), (40, 48))         # far away: out of scope
    row, pt = n3.benchmark_gonad(pred, truth, vox, coverage="protocol", min_truth_voxels=10)
    assert row["coverage_rule"].startswith("protocol") and row["n_truth"] == 2
    assert row["n_truth_fragments_removed"] == 1 and row["truth_fragments_removed"] == "3"
    # region C reaches 1.0 um = 10 voxels in y beyond the masks (rows 15 to 24): prediction 13 has rows
    # 20 to 25, five of its six rows inside C, so it is in scope; prediction 14 is not
    assert row["n_pred_ignored_outside_coverage"] == 1 and row["n_pred_counted"] == 4
    assert row["tp_iou50"] == 2 and row["fp_iou50"] == 2 and row["fn_iou50"] == 0
    # 20 percent rule: truth 1 is covered by 10 (70 percent) and 11 (30 percent): a split
    assert row["n_splits"] == 1 and row["n_merges"] == 0 and row["split_rule"].startswith("in_scope_cover")
    assert pt.set_index("truth_label").loc[1, "is_split"] and pt.set_index("truth_label").loc[1, "n_pred_split_rule"] == 2
    # the majority rule of the box modes says otherwise for prediction 11 (it lies wholly in truth 1)
    row_auto, _ = n3.benchmark_gonad(pred, truth, vox, coverage="bbox")
    assert row_auto["n_splits"] == 1 and row_auto["n_truth_fragments_removed"] == 0
    # a merge: one prediction covering 25 percent of two nuclei
    pred2 = np.zeros_like(truth)
    _cube(pred2, 5, (2, 6), (5, 15), (12, 18))
    row2, _ = n3.benchmark_gonad(pred2, truth, vox, coverage="protocol", min_truth_voxels=10)
    assert row2["n_merges"] == 1
    # the z-extent error in slices
    pred3 = truth.copy()
    pred3[pred3 == 3] = 0
    pred3[1][truth[2] == 1] = 1                          # truth 1 one plane taller
    row3, pt3 = n3.benchmark_gonad(pred3, truth, vox, coverage="protocol", min_truth_voxels=10)
    assert abs(row3["z_extent_error_mean_slices"] - 0.5) < 1e-9 and abs(row3["z_extent_error_mean_um"] - 0.1) < 1e-9
    assert pt3["z_extent_error_slices"].tolist() == [1.0, 0.0]
    with pytest.raises(ValueError, match="coverage must be"):
        n3.benchmark_gonad(pred, truth, vox, coverage="nope")


def test_protocol_scope_is_half_the_volume():
    vox = (0.2, 0.1, 0.1)
    # (a 400-voxel truth object: the fragment floor is lowered from the protocol's 1000)
    truth = np.zeros((4, 40, 80), np.int32)
    _cube(truth, 1, (1, 3), (10, 20), (10, 20))
    pred = np.zeros_like(truth)
    # beside the truth in x, in its rows and planes: C reaches x = 29 (1.0 um = 10 voxels beyond x = 19),
    # so 10 of the 20 columns are inside: exactly half, in scope
    _cube(pred, 1, (1, 3), (10, 20), (20, 40))
    row, _ = n3.benchmark_gonad(pred, truth, vox, coverage="protocol", min_truth_voxels=10)
    assert row["n_pred_counted"] == 1 and row["fp_iou50"] == 1
    pred[pred == 1] = 0
    _cube(pred, 1, (1, 3), (10, 20), (21, 41))           # 9 of 20 columns inside: out of scope
    row, _ = n3.benchmark_gonad(pred, truth, vox, coverage="protocol", min_truth_voxels=10)
    assert row["n_pred_counted"] == 0 and row["fp_iou50"] == 0 and row["fn_iou50"] == 1


def test_benchmark_script_scores_a_cellpose_model_file(tiny_bench, monkeypatch):
    """--cellpose-model runs the pipeline's Cellpose path (stubbed here) on the cached DAPI and scores it."""
    from scipy import ndimage as ndi

    from germquant.segment import nuclei as N

    calls = []

    def fake_cellpose(dna, spacing, model_name, diameter_um, device_pref=None):
        calls.append(str(model_name))
        return ndi.label(np.asarray(dna) > 500)[0].astype(np.int32)

    monkeypatch.setattr(N, "_cellpose", fake_cellpose)
    b = tiny_bench
    model = b["root"] / "cp_model"
    model.write_bytes(b"a stand-in for a Cellpose model file")
    spec = importlib.util.spec_from_file_location("benchmark_nuclei3d", ROOT / "scripts" / "benchmark_nuclei3d.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    out = b["root"] / "bench_cp"
    rc = mod.main(["--manifest", str(b["manifest"]), "--out", str(out), "--coverage", "protocol", "--split", "val",
                   "--cellpose-model", str(model), "--cache-root", str(b["root"] / "cache"), "--device", "cpu",
                   "--n-boot", "20"])
    assert rc == 0 and calls == [str(model)]
    row = pd.read_csv(out / "per_gonad.csv").iloc[0]
    assert row["pred_source"] == "cellpose model" and row["segmentation_method"] == "cellpose"
    assert row["tp_iou50"] == 4 and row["fp_iou50"] == 0 and row["f1_iou50"] == 1.0
    pooled = json.loads((out / "pooled.json").read_text(encoding="utf-8"))
    assert pooled["settings"]["model"]["kind"] == "cellpose"
    assert (out / "labels" / "synth_val_c__labels.tif").is_file()


def test_export_imaris_training_pairs(tiny_bench):
    """The Cellpose competitor's slice export: train gonads only, the truth masks as labels, one series."""
    import tifffile

    from germquant.segment.training import load_pairs, split_pairs

    b = tiny_bench
    spec = importlib.util.spec_from_file_location("export_pairs", ROOT / "scripts" / "export_imaris_training_pairs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    out = b["root"] / "pairs"
    rc = mod.main(["--out", str(out), "--split-file", str(b["split"]), "--manifest", str(b["manifest"]),
                   "--cache-root", str(b["root"] / "cache"), "--slices-per-image", "2", "--min-labels", "1"])
    assert rc == 0
    man = pd.read_csv(out / "manifest.csv")
    assert set(man["image_id"]) == {"synth_noHS_a", "synth_HS_b"} and (man["series_id"] == "train").all()
    assert set(man["orientation"]) == {"xy", "xz", "yz"} and len(man) == 12
    exp = json.loads((out / "export.json").read_text(encoding="utf-8"))
    assert exp["test_never_opened"] == ["synth_test_d"] and exp["mask_untraced"] is True
    pairs = load_pairs(out)
    train, val, rec = split_pairs(pairs, by="series")
    assert len(val) == 0 and len(train) == len(man)            # a single group: everything trains
    for r in man.itertuples():
        lab = tifffile.imread(out / r.labels)
        im = tifffile.imread(out / r.image)
        assert lab.shape == im.shape == (r.height, r.width) and im.dtype == np.uint16
        assert set(np.unique(lab).tolist()) <= {0, 1, 2, 3, 4} and lab.max() > 0
        if r.orientation != "xy":
            assert abs(r.z_zoom - SP[0] / SP[1]) < 1e-9
