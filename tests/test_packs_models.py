"""Pack models: the blob store, fetch-models, resolution and the pin check (design 3.5, P1 acceptance 4 and 5).

Small fakes stand in for the 1.2 GB weights and for the model server: the test pack's stand-in model
tinynet@f5055f8a (64 bytes) and a local HTTP server with Range support. fetch-models adopts a Cellpose-cache
file and a --from-file copy after size and sha256 checks, downloads and resumes, discards a mismatch and does
nothing on a second call; the blob lands in models_dir()/blobs/<sha256>/ and resolve_model_path returns it; a
`pack check` in a new process trusts a blob whose size and mtime match its record and re-hashes after a change.
A local registry entry with a pack model's name but other weights is refused before segmentation; names of the
pack form are reserved; a declared but missing blob raises the fetch-models message under cellpose and auto,
while classical passes resolution and records model_sha256 "unused".
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import time

import numpy as np
import pytest
import yaml
from packs_helpers import (  # noqa: F401  (pack_env is a fixture)
    PYTHON,
    TEST_PACK,
    TINYNET,
    TINYNET_BYTES,
    TINYNET_SHA,
    copy_pack,
    pack_env,
    patch_nd2_reader,
    small_stack,
    subprocess_env,
    use_fixture_packs,
)
from test_datafetch import serve

from germquant import models_registry as MR
from germquant import packs, pipeline
from germquant.config import load_config
from germquant.segment import backends as SB
from germquant.segment import nuclei as N

NAME = "20261011_field_001.nd2"


def _blob(env):
    return env.blobs / TINYNET_SHA / "tinynet"


def _with_model_cfg(method=None):
    cfg = load_config(TEST_PACK / "config" / "profiles" / "with_model.yaml")
    if method:
        cfg._data["segmentation"]["nuclei"]["method"] = method
    return cfg


def test_fetch_adopts_the_cellpose_cache_and_is_idempotent(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    (pack_env.cache / "tinynet").write_bytes(TINYNET_BYTES)
    res = packs.fetch_pack_models("test_pack", "with_model")
    assert [(r["name"], r["status"]) for r in res] == [(TINYNET, "adopted")]
    blob = _blob(pack_env)
    assert blob.read_bytes() == TINYNET_BYTES
    rec = json.loads((pack_env.blobs / TINYNET_SHA / f"{TINYNET_SHA}.json").read_text(encoding="utf-8"))
    assert rec["sha256"] == TINYNET_SHA and rec["size"] == 64 and "cellpose_cache" in rec["source"]
    assert rec["packs"] == ["test_pack@0.1.0"] and rec["verified_at"]
    mtime = blob.stat().st_mtime_ns
    again = packs.fetch_pack_models("test_pack", "with_model")
    assert again[0]["status"] == "present" and blob.stat().st_mtime_ns == mtime
    # resolution: the pinned name resolves to the verified blob, through the pipeline's own lookup
    assert pipeline.resolve_model_path(_with_model_cfg()) == blob
    assert MR.resolve_model_name(TINYNET) == blob
    assert MR.resolve_model_name(TINYNET, expect_sha256=TINYNET_SHA) == blob
    with pytest.raises(ValueError, match="expected"):
        MR.resolve_model_name(TINYNET, expect_sha256="0" * 64)


def test_fetch_from_file_and_mismatch_is_never_registered(pack_env, monkeypatch, tmp_path):  # noqa: F811
    use_fixture_packs(monkeypatch)
    wrong = TINYNET_BYTES.replace(b"tinynet", b"tinyNET")             # the same size, other bytes
    (pack_env.cache / "tinynet").write_bytes(wrong)
    bad_file = tmp_path / "bad.bin"
    bad_file.write_bytes(b"x" * 64)
    res = packs.fetch_pack_models("test_pack", "with_model", {TINYNET: bad_file}, download=False)
    assert res[0]["status"] == "missing" and not _blob(pack_env).exists()
    assert any("sha256" in t for t in res[0]["tried"])               # the cache copy failed its check
    good = tmp_path / "tinynet.bin"
    good.write_bytes(TINYNET_BYTES)
    res = packs.fetch_pack_models("test_pack", "with_model", {TINYNET: good})
    assert res[0]["status"] == "from_file" and _blob(pack_env).read_bytes() == TINYNET_BYTES
    with pytest.raises(ValueError, match="does not declare"):
        packs.fetch_pack_models("test_pack", None, {"other@12345678": good})


def _url_pack(tmp_path, url):
    def edit(d):
        d["models"][0]["sources"] = [{"type": "url", "url": url}]
    return copy_pack(TEST_PACK, tmp_path / "url_packs" / "test_pack" / "0.1.0", edit=edit)


def test_fetch_downloads_resumes_and_discards(pack_env, monkeypatch, tmp_path):  # noqa: F811
    with serve(payload=TINYNET_BYTES, cut_first_at=20) as (url, h):
        _url_pack(tmp_path, url)
        monkeypatch.setenv("GERMQUANT_PACKS", str(tmp_path / "url_packs"))
        packs.refresh()
        first = packs.fetch_pack_models("test_pack", "with_model")
        assert first[0]["status"] == "missing"                          # cut short: the partial file stays
        part = _blob(pack_env).with_name("tinynet.part")
        assert part.is_file() and part.stat().st_size == 20
        second = packs.fetch_pack_models("test_pack", "with_model")
        assert second[0]["status"] == "downloaded" and h.log[-1][1] == "bytes=20-"
    assert _blob(pack_env).read_bytes() == TINYNET_BYTES and not part.exists()
    # a server with other bytes: discarded, never registered
    _blob(pack_env).unlink()
    (pack_env.blobs / TINYNET_SHA / f"{TINYNET_SHA}.json").unlink()
    with serve(payload=b"y" * 64) as (url, h):
        shutil.rmtree(tmp_path / "url_packs")
        _url_pack(tmp_path, url)
        packs.refresh()
        res = packs.fetch_pack_models("test_pack", "with_model")
    assert res[0]["status"] == "missing" and not _blob(pack_env).exists() and not part.exists()


def _check_json(env_vars) -> dict:
    out = subprocess.run([PYTHON, "-m", "germquant.packs", "check", "test_pack", "--recipe", "with_model", "--json"],
                         capture_output=True, text=True, env=subprocess_env(**env_vars), timeout=300)
    return json.loads(out.stdout)


def test_check_in_a_new_process_trusts_the_record_and_rehashes_after_a_change(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    (pack_env.cache / "tinynet").write_bytes(TINYNET_BYTES)
    packs.fetch_pack_models("test_pack", "with_model")
    env = {"GERMQUANT_PACKS": os.environ["GERMQUANT_PACKS"], "GERMQUANT_MODELS": os.environ["GERMQUANT_MODELS"],
           "LOCALAPPDATA": os.environ["LOCALAPPDATA"], "XDG_DATA_HOME": os.environ["XDG_DATA_HOME"]}
    st = _check_json(env)
    assert st["models"][0]["state"] == "trusted", st                    # size and mtime match: not re-hashed
    blob = _blob(pack_env)
    later = time.time() + 120
    os.utime(blob, (later, later))
    st = _check_json(env)
    assert st["models"][0]["state"] == "verified", st                   # mtime moved: re-hashed, record refreshed
    assert _check_json(env)["models"][0]["state"] == "trusted"
    blob.write_bytes(TINYNET_BYTES.replace(b"tinynet", b"tinyNET"))
    st = _check_json(env)
    assert st["models"][0]["state"] == "mismatch" and not st["usable"]
    assert any("fetch-models test_pack" in b["fix"] for b in st["blockers"])


def test_missing_blob_raises_the_fetch_models_message_except_under_classical(pack_env, monkeypatch, tmp_path):  # noqa: F811
    use_fixture_packs(monkeypatch)
    for method in ("cellpose", "auto"):
        with pytest.raises(FileNotFoundError, match="germquant pack fetch-models test_pack"):
            pipeline.resolve_model_path(_with_model_cfg(method))
    cfg = _with_model_cfg("classical")
    assert pipeline.resolve_model_path(cfg) == TINYNET                  # returned unresolved
    patch_nd2_reader(monkeypatch, small_stack())
    res = pipeline.process_image(NAME, cfg, tmp_path / "out", prov=None)
    man = json.loads((tmp_path / "out" / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["model_sha256"] == "unused" and res["n_nuclei"] == 5
    # a name no installed pack declares
    cfg2 = _with_model_cfg("cellpose")
    cfg2._data["segmentation"]["nuclei"]["cellpose_model"] = "ghost@0123abcd"
    with pytest.raises(FileNotFoundError, match="no installed or built-in pack declares it"):
        pipeline.resolve_model_path(cfg2)


def test_local_entry_with_a_pack_name_but_other_weights_is_refused_before_segmentation(
        pack_env, monkeypatch, tmp_path):  # noqa: F811
    use_fixture_packs(monkeypatch)
    other = tmp_path / "other_weights"
    other.write_bytes(b"z" * 64)
    reg = MR.local_registry_path()
    reg.parent.mkdir(parents=True, exist_ok=True)
    reg.write_text(json.dumps({"format": 1, "models": {TINYNET: {"path": str(other), "sha256": "z", "size": 64}}}),
                   encoding="utf-8")                                    # an entry written before the name was reserved
    calls = []
    monkeypatch.setattr(SB, "_REGISTRY", dict(SB._REGISTRY))
    monkeypatch.setattr(SB, "_SOURCES", dict(SB._SOURCES))     # restored afterwards, as both dicts must be
    SB.register_backend("cellpose", lambda *a, **k: calls.append(1) or np.zeros(a[0].shape, np.int32), replace=True)
    patch_nd2_reader(monkeypatch, small_stack())
    with pytest.raises(MR.ModelPinError, match="does not start with its pin f5055f8a"):
        pipeline.process_image(NAME, _with_model_cfg(), tmp_path / "out", prov=None)
    assert calls == []                                                  # segmentation never started
    assert not (tmp_path / "out" / "run_manifest.json").exists()


def test_pin_check_compares_with_the_spec_digest(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    MR.check_model_pin(TINYNET, TINYNET_SHA)
    MR.check_model_pin("models/models/germline_nuclei_combined", "anything")      # not a pack name: no check
    with pytest.raises(MR.ModelPinError, match="the pack spec pins"):
        MR.check_model_pin(TINYNET, "f5055f8a" + "0" * 56)                     # right prefix, wrong digest
    with pytest.raises(MR.ModelPinError, match="does not start"):
        MR.check_model_pin(TINYNET, "absent")
    # a pinned name no pack declares: 8 hex digits alone pin nothing, so it is refused
    with pytest.raises(MR.ModelPinError, match="no installed or built-in pack declares it"):
        MR.check_model_pin("ghost@0123abcd", "0123abcd" + "f" * 56)


def test_pack_names_are_reserved(pack_env, tmp_path):  # noqa: F811
    f = tmp_path / "m"
    f.write_bytes(b"w" * 10)
    with pytest.raises(ValueError, match="reserved for pack models"):
        MR.register_local_model(f, name="mymodel@0123abcd")
    rec = MR.register_local_model(f, name="lab@v2")                  # an '@' alone is not the pack form
    assert rec["name"] == "lab@v2" and MR.resolve_registered("lab@v2") == f.resolve()
    assert MR.is_pack_model_name("cpsam_v2@0f1cc3f7") and not MR.is_pack_model_name("cpsam_v2")
    assert not MR.is_pack_model_name("models/x@0123abcd") and not MR.is_pack_model_name("lab@v2")
    for builtin in MR.CELLPOSE_BUILTIN_NAMES + tuple(MR.MODELS):
        assert not MR.is_pack_model_name(builtin)


def test_fetch_discards_a_wrong_blob_already_in_the_store(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    blob = _blob(pack_env)
    blob.parent.mkdir(parents=True)
    blob.write_bytes(TINYNET_BYTES.replace(b"tinynet", b"tinyNET"))            # the right size, other bytes
    res = packs.fetch_pack_models("test_pack", "with_model", download=False)
    assert res[0]["status"] == "missing" and not blob.exists()                 # discarded, never kept
    assert any(t.startswith("store: sha256") for t in res[0]["tried"])


def test_a_candidate_of_the_wrong_size_is_never_hashed(pack_env, monkeypatch, tmp_path):  # noqa: F811
    from germquant import datafetch

    use_fixture_packs(monkeypatch)
    (pack_env.cache / "tinynet").write_bytes(TINYNET_BYTES[:-1])
    longer = tmp_path / "longer.bin"
    longer.write_bytes(TINYNET_BYTES + b"x")
    hashed = []
    real = datafetch.sha256_file
    monkeypatch.setattr(datafetch, "sha256_file", lambda p, *a, **k: hashed.append(str(p)) or real(p, *a, **k))
    res = packs.fetch_pack_models("test_pack", "with_model", {TINYNET: longer}, download=False)
    assert res[0]["status"] == "missing" and hashed == []                      # the size failed first
    assert any("63 bytes, expected 64" in t for t in res[0]["tried"])
    assert any("65 bytes, expected 64" in t for t in res[0]["tried"])


def test_check_blob_verifies_when_the_store_cannot_keep_a_record(pack_env, monkeypatch, caplog):  # noqa: F811
    from germquant.packs import blobs as BL

    use_fixture_packs(monkeypatch)
    # a short store root: a folder named like the record, two sha256 deep, passes Windows' 248-character
    # limit for folders only from a short path
    with tempfile.TemporaryDirectory(prefix="gqm") as short:
        monkeypatch.setenv("GERMQUANT_MODELS", short)
        blob = BL.blobs_dir() / TINYNET_SHA / "tinynet"
        blob.parent.mkdir(parents=True)
        blob.write_bytes(TINYNET_BYTES)
        rec = BL.record_path(TINYNET_SHA)
        rec.mkdir()                                   # the record cannot be written (as on a read-only store)
        spec = packs.find_pack("test_pack").model(TINYNET)
        assert BL.check_blob(spec)["state"] == "verified"
        assert "could not save the verification record" in caplog.text
        assert not rec.with_name(rec.name + ".part").exists()                   # no stale partial record
        st = packs.pack_status("test_pack", recipe="with_model")
        assert st["models"][0]["state"] == "verified" and not any("tinynet" in b["what"] for b in st["blockers"])
        assert packs.fetch_pack_models("test_pack", "with_model", download=False)[0]["status"] == "present"
        monkeypatch.setenv("GERMQUANT_MODELS", str(pack_env.models))


def test_cellpose_backend_first_call_is_the_old_line(monkeypatch):
    """Rule R3: on success cellpose_backend makes exactly the old call, once, and returns its labels as they
    are; nothing else runs (no device lookup, no record). Only a CPU retry records "cpu"."""
    from germquant import device

    calls = []
    labels = np.ones((2, 4, 4), np.int32)

    def fake_cellpose(dna, spacing, model, diam, device_pref=None):
        calls.append((model, diam, device_pref))
        if device_pref != "cpu" and fail_first:
            raise RuntimeError("CUDA out of memory")
        return labels

    def no_lookup(*a, **k):
        raise AssertionError("the success path looked up the device")

    monkeypatch.setattr(N, "_cellpose", fake_cellpose)
    monkeypatch.setattr(device, "select_device", no_lookup)
    device.record_segment_device(None)
    fail_first = False
    dna = np.zeros((2, 4, 4), np.float32)
    out = SB.cellpose_backend(dna, (1.0, 0.2, 0.2), model="models/m", diameter_um=3.0, device="cuda", cfg_section={})
    assert out is labels and calls == [("models/m", 3.0, "cuda")] and device.segment_device_used() is None
    monkeypatch.undo()
    monkeypatch.setattr(N, "_cellpose", fake_cellpose)
    monkeypatch.setattr(device, "available", lambda: {"cuda": True, "mps": False})
    monkeypatch.delenv("GERMQUANT_DEVICE", raising=False)
    calls.clear()
    fail_first = True
    out = SB.cellpose_backend(dna, (1.0, 0.2, 0.2), model="models/m", diameter_um=3.0, device="cuda", cfg_section={})
    assert out is labels and [c[2] for c in calls] == ["cuda", "cpu"] and device.segment_device_used() == "cpu"
    device.record_segment_device(None)


def test_train_base_refuses_a_local_entry_with_a_pack_name_but_other_weights(pack_env, monkeypatch, tmp_path):  # noqa: F811
    from germquant.segment import training as TR

    use_fixture_packs(monkeypatch)
    reg = MR.local_registry_path()
    reg.parent.mkdir(parents=True, exist_ok=True)
    other, good = tmp_path / "other_weights", tmp_path / "good_weights"
    other.write_bytes(b"z" * 64)
    good.write_bytes(TINYNET_BYTES)
    for path in (other, good):                        # entries written before the name was reserved
        reg.write_text(json.dumps({"format": 1, "models": {TINYNET: {"path": str(path), "size": 64}}}),
                       encoding="utf-8")
        if path is other:
            with pytest.raises(MR.ModelPinError, match="does not start with its pin f5055f8a"):
                TR._resolve_base(TINYNET)
        else:
            assert TR._resolve_base(TINYNET) == {"name": TINYNET, "kind": "registered", "path": str(good),
                                                 "sha256": TINYNET_SHA}


def test_run_records_the_blob_digest_and_the_segment_device(pack_env, monkeypatch, tmp_path):  # noqa: F811
    """The real cellpose_backend path with nuclei._cellpose faked: the manifest's model_sha256 is the blob's
    digest, and the pack run records the device Cellpose ran on (CPU here, CUDA hidden)."""
    from germquant import api, device

    use_fixture_packs(monkeypatch)
    (pack_env.cache / "tinynet").write_bytes(TINYNET_BYTES)
    packs.fetch_pack_models("test_pack", "with_model")
    seen = []

    def fake_cellpose(dna, spacing, model, diam, device_pref=None):
        seen.append(str(model))
        lab, _ = N.segment_nuclei(dna, spacing, method="classical", diameter_um=3.0, min_volume_um3=2.0)
        return lab

    monkeypatch.setattr(N, "_cellpose", fake_cellpose)
    monkeypatch.setattr(device, "available", lambda: {"cuda": False, "mps": False})
    patch_nd2_reader(monkeypatch, small_stack())
    res = api.run(NAME, pack="test_pack", recipe="with_model", out=tmp_path / "out")
    assert seen == [str(_blob(pack_env))] and res["n_nuclei"] == 5
    out = tmp_path / "out"
    man = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["model_sha256"] == TINYNET_SHA and man["pack"]["recipe"] == "with_model"
    rec = json.loads((out / f"{res['image_id']}__stages.json").read_text(encoding="utf-8"))
    assert rec["post_run"] == {"segment_device_used": "cpu"}
    assert rec["runtime"]["device_selected"] == "cpu" and rec["model_sha256"] == TINYNET_SHA


def test_train_base_resolves_a_pack_model(pack_env, monkeypatch, tmp_path):  # noqa: F811
    from germquant.segment import training as TR

    use_fixture_packs(monkeypatch)
    with pytest.raises(FileNotFoundError, match="fetch-models"):
        TR._resolve_base(TINYNET)
    good = tmp_path / "t.bin"
    good.write_bytes(TINYNET_BYTES)
    packs.fetch_pack_models("test_pack", "with_model", {TINYNET: good})
    base = TR._resolve_base(TINYNET)
    assert base == {"name": TINYNET, "kind": "pack", "path": str(_blob(pack_env)), "sha256": TINYNET_SHA}
    assert TR._resolve_base("cpsam")["kind"] == "cellpose_builtin"           # unchanged


def test_manifest_spec_names_the_model_and_card(pack_env):  # noqa: F811
    data = yaml.safe_load((TEST_PACK / "pack.yaml").read_text(encoding="ascii"))
    spec = data["models"][0]
    assert spec["name"] == TINYNET and spec["sha256"] == TINYNET_SHA and spec["size"] == len(TINYNET_BYTES)
    card = json.loads((TEST_PACK / spec["card"]).read_text(encoding="utf-8"))
    assert (card["stand_in_content"] + "\n").encode("ascii") == TINYNET_BYTES
