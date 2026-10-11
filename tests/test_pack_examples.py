"""`germquant pack fetch-example`, germquant.packs.fetch_example and api.fetch_example (design 3.3, 3.9, 3.13 and 6.9;
package P3): a pack's real example is fetched from a local stand-in server (plain files, and one member of a zip
read with HTTP Range requests), checked by size and sha256, assembled into one OME-TIFF and checked by
pixels_sha256; a second call downloads nothing; a wrong digest is refused and nothing wrong is kept; the
torch-free entry `python -m germquant.packs fetch-example` loads no torch, CuPy, SpotMAX or Cellpose; malformed
example specs make a pack unusable. Nothing leaves the machine: the specs name an https address that the tests
map to the local server."""
from __future__ import annotations

import hashlib
import http.server
import io
import json
import os
import re
import subprocess
import threading
import zipfile
from contextlib import contextmanager

import numpy as np
import pytest
import tifffile
import yaml
from packs_helpers import HEAVY, PYTHON, ROOT, TEST_PACK, copy_pack, pack_env, subprocess_env  # noqa: F401

from germquant import api, cli_packs, datafetch, packs
from germquant.assemble import build, pixels_sha256, read_ome_tiff
from germquant.packs import examples as E

BASE = "https://data.example.invalid/files/"
SP = [0.5, 0.2, 0.2]
REAL_SOURCE_URL = E.source_url            # the fixture maps the https address to the local server


class _Files(http.server.BaseHTTPRequestHandler):
    files: dict = {}
    log: list = []

    def log_message(self, *a):
        pass

    def do_GET(self):  # noqa: N802 - the http.server API
        cls = type(self)
        path = self.path.lstrip("/")
        cls.log.append((path, self.headers.get("Range")))
        data = cls.files.get(path)
        if data is None:
            self.send_response(404)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        rng = self.headers.get("Range")
        start, end = 0, len(data) - 1
        if rng:
            a, b = re.fullmatch(r"bytes=(\d*)-(\d*)", rng).groups()
            start, end = (len(data) - int(b), len(data) - 1) if a == "" else (int(a), int(b) if b else len(data) - 1)
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
        else:
            self.send_response(200)
        body = data[start:end + 1]
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


@contextmanager
def serve(files: dict):
    handler = type("H", (_Files,), {"files": files, "log": []})
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        yield f"http://127.0.0.1:{srv.server_address[1]}/", handler
    finally:
        srv.shutdown()
        srv.server_close()


def _tif_bytes(arr, **kw) -> bytes:
    buf = io.BytesIO()
    tifffile.imwrite(buf, arr, **kw)
    return buf.getvalue()


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


@pytest.fixture
def remote(tmp_path, monkeypatch, pack_env):  # noqa: F811
    """A copy of the test pack (renamed fetch_pack) with two real examples served locally: two per-channel TIFFs, and
    one member of a zip holding a multi-channel TIFF with channels as samples (DeepFoci's layout)."""
    rng = np.random.default_rng(7)
    dna = rng.integers(100, 3000, (6, 32, 40), dtype=np.uint16)
    foci = rng.integers(100, 3000, (6, 32, 40), dtype=np.uint16)
    rgb = rng.integers(100, 3000, (6, 32, 40, 3), dtype=np.uint16)
    files = {"set/x_dna.tif": _tif_bytes(dna, imagej=True), "set/x_foci.tif": _tif_bytes(foci)}
    zbuf = io.BytesIO()
    with zipfile.ZipFile(zbuf, "w") as zf:
        zf.writestr("valid/readme.txt", b"stand-in archive")
        zf.writestr(zipfile.ZipInfo("valid/data_001.tif"), _tif_bytes(rgb, photometric="rgb"))   # stored
    files["record/archive.zip/content"] = zbuf.getvalue()
    member = _tif_bytes(rgb, photometric="rgb")
    root = tmp_path / "packs_root"

    def spec_planes(**over):
        tmp_inputs = tmp_path / "inputs"
        tmp_inputs.mkdir(exist_ok=True)
        for rel in ("set/x_dna.tif", "set/x_foci.tif"):
            (tmp_inputs / rel.split("/")[1]).write_bytes(files[rel])
        data, _ = build([tmp_inputs / "x_dna.tif", tmp_inputs / "x_foci.tif"], ["DAPI", "foci"], axes="ZYX",
                        crop={"y": [2, 30]})
        spec = {"id": "planes", "base_url": BASE,
                "sources": [{"path": rel, "size": len(files[rel]), "sha256": _sha(files[rel])}
                            for rel in ("set/x_dna.tif", "set/x_foci.tif")],
                "assemble": {"axes": "ZYX", "channels": ["DAPI", "foci"], "spacing_zyx_um": SP, "crop": {"y": [2, 30]}},
                "pixels_sha256": pixels_sha256(data), "licence": "CC0", "attribution": "stand-in data, test only",
                "expect": {"n_objects": [1, 9], "spots_per_object": [0.5, 3.0]}}
        spec.update(over)
        return spec

    def spec_member(**over):
        (tmp_path / "inputs").mkdir(exist_ok=True)
        (tmp_path / "inputs" / "data_001.tif").write_bytes(member)
        data, _ = build([tmp_path / "inputs" / "data_001.tif"], ["53BP1", "foci", "DAPI"], axes="ZYXC")
        zinfo = zipfile.ZipFile(io.BytesIO(files["record/archive.zip/content"])).getinfo("valid/data_001.tif")
        spec = {"id": "member", "base_url": BASE,
                "sources": [{"path": "record/archive.zip/content", "member": "valid/data_001.tif",
                             "zip_size": len(files["record/archive.zip/content"]), "crc32": f"{zinfo.CRC:08x}",
                             "size": len(member), "sha256": _sha(member)}],
                "assemble": {"axes": "ZYXC", "channels": ["53BP1", "foci", "DAPI"], "spacing_zyx_um": SP},
                "pixels_sha256": pixels_sha256(data), "licence": "CC BY 4.0", "attribution": "stand-in, test only"}
        spec.update(over)
        return spec

    def make(*specs):
        if (root / "fetch_pack").exists():
            import shutil

            shutil.rmtree(root / "fetch_pack")

        def edit(d):
            d["id"] = "fetch_pack"
            d["examples"] = [{"id": "synthetic", "image": "example/test_pack_synthetic.ome.tif"}, *specs]
        copy_pack(TEST_PACK, root / "fetch_pack", edit=edit)
        monkeypatch.setenv("GERMQUANT_PACKS", str(root))
        packs.refresh()
        return packs.find_pack("fetch_pack")

    with serve(files) as (base, handler):
        monkeypatch.setattr(E, "source_url", lambda spec, src: base + str(src["path"]))
        yield type("R", (), {"make": staticmethod(make), "planes": staticmethod(spec_planes),
                             "member": staticmethod(spec_member), "log": handler.log, "base": base,
                             "files": files, "out": tmp_path / "examples", "dna": dna, "foci": foci, "rgb": rgb})


def test_fetch_assemble_record_and_a_second_call_downloads_nothing(remote):
    pack = remote.make(remote.planes(), remote.member())
    rec = packs.fetch_example(pack, "planes", remote.out)
    assert rec["status"] == "fetched" and rec["kind"] == "remote" and rec["recipe"] == "basic"
    img = read_ome_tiff(rec["image"])
    assert img["channel_names"] == ["DAPI", "foci"] and img["spacing_zyx_um"] == pytest.approx(SP)
    assert np.array_equal(img["data"], np.stack([remote.dna, remote.foci])[:, :, 2:30])
    assert img["pixels_sha256"] == rec["pixels_sha256"]
    record = json.loads((remote.out / "fetch_pack" / "planes" / "planes.json").read_text(encoding="ascii"))
    assert record["licence"] == "CC0" and record["attribution"] == "stand-in data, test only"
    assert [s["status"] for s in record["sources"]] == ["downloaded", "downloaded"]
    assert record["shape_czyx"] == [2, 6, 28, 40] and record["crop"] == {"y": [2, 30]}
    n_requests = len(remote.log)
    again = packs.fetch_example(pack, "planes", remote.out)
    assert again["status"] == "present" and len(remote.log) == n_requests          # nothing downloaded
    for f in (remote.out / "fetch_pack" / "planes" / "source").iterdir():          # sources deleted to save space
        f.unlink()
    assert packs.fetch_example(pack, "planes", remote.out)["status"] == "present" and len(remote.log) == n_requests
    (remote.out / "fetch_pack" / "planes" / "planes.json").unlink()                # an image without its record
    third = packs.fetch_example(pack, "planes", remote.out)
    assert third["status"] == "present" and [s["status"] for s in third["sources"]] == ["not kept", "not kept"]
    rewritten = json.loads((remote.out / "fetch_pack" / "planes" / "planes.json").read_text(encoding="ascii"))
    assert rewritten["licence"] == "CC0" and rewritten["shape_czyx"] == [2, 6, 28, 40]
    assert rewritten["channel_names"] == ["DAPI", "foci"] and len(remote.log) == n_requests
    # one member of a zip, read with Range requests only
    rec2 = packs.fetch_example(pack, "member", remote.out)
    member_requests = [r for p, r in remote.log[n_requests:] if p == "record/archive.zip/content"]
    assert member_requests and all(r and r.startswith("bytes=") for r in member_requests)
    img2 = read_ome_tiff(rec2["image"])
    assert img2["channel_names"] == ["53BP1", "foci", "DAPI"]
    assert np.array_equal(img2["data"], np.moveaxis(remote.rgb, -1, 0))
    assert api.fetch_example is packs.fetch_example and "fetch_example" in api.__all__


def test_the_default_example_is_the_first_real_one_and_a_local_one_is_in_the_pack(remote):
    pack = remote.make(remote.planes())
    assert packs.find_example(pack)["id"] == "planes"
    rec = packs.fetch_example(pack, "synthetic", remote.out)
    assert rec["status"] == "local" and rec["image"] == str(pack.folder / "example" / "test_pack_synthetic.ome.tif")
    with pytest.raises(FileNotFoundError, match="no example 'nope'.*synthetic, planes"):
        packs.fetch_example(pack, "nope", remote.out)


def test_wrong_pins_are_refused_and_nothing_wrong_is_kept(remote):
    pack = remote.make(remote.planes(pixels_sha256="0" * 64))
    with pytest.raises(packs.ExampleError, match="pixels have sha256"):
        packs.fetch_example(pack, "planes", remote.out)
    folder = remote.out / "fetch_pack" / "planes"
    assert not (folder / "planes.ome.tif").exists() and not (folder / "planes.json").exists()
    assert (folder / "source" / "x_dna.tif").is_file()                     # the sources matched, so they stay
    bad = remote.planes()
    bad["sources"][1]["sha256"] = "f" * 64
    pack = remote.make(bad)
    with pytest.raises(datafetch.DownloadError, match="sha256"):
        packs.fetch_example(pack, "planes", remote.out / "b")
    assert not (remote.out / "b" / "fetch_pack" / "planes" / "source" / "x_foci.tif").exists()
    with pytest.raises(packs.ExampleError, match="no-download"):
        packs.fetch_example(remote.make(remote.planes()), "planes", remote.out / "c", download=False)


def test_fetch_example_command_and_its_torch_free_entry(remote, capsys):
    remote.make(remote.planes())
    assert cli_packs.main(["fetch-example", "fetch_pack", "--example", "planes", "--out", str(remote.out)]) == 0
    out = capsys.readouterr().out
    assert "fetched" in out and "--pack fetch_pack@0.1.0" in out and "n_objects in [1, 9]" in out
    assert cli_packs.main(["fetch-example", "fetch_pack", "--out", str(remote.out), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "present"
    assert cli_packs.main(["fetch-example", "fetch_pack", "--example", "x", "--out", str(remote.out)]) == 1
    # a fresh interpreter: the light entry fetches and assembles without torch, CuPy, SpotMAX or Cellpose
    code = (
        "import json, sys\n"
        "from germquant.packs import examples as E\n"
        f"E.source_url = lambda spec, src: {remote.base!r} + str(src['path'])\n"
        "from germquant.cli_packs import main\n"
        f"rc = main(['fetch-example', 'fetch_pack', '--example', 'planes', '--out', {str(remote.out / 'fresh')!r}])\n"
        f"print(json.dumps({{'rc': rc, 'heavy': [m for m in {list(HEAVY)!r} if m in sys.modules]}}))\n")
    env = subprocess_env(GERMQUANT_PACKS=os.environ["GERMQUANT_PACKS"], GERMQUANT_MODELS=os.environ["GERMQUANT_MODELS"])
    r = subprocess.run([PYTHON, "-c", code], capture_output=True, text=True, env=env, timeout=300)
    assert r.returncode == 0, r.stderr
    res = json.loads(r.stdout.strip().splitlines()[-1])
    assert res == {"rc": 0, "heavy": []}, r.stdout
    assert (remote.out / "fresh" / "fetch_pack" / "planes" / "planes.ome.tif").is_file()


def test_a_loopback_address_serves_tests_and_lint_refuses_it(remote, monkeypatch):
    """Plain http to this machine is accepted, so a stand-in server works across processes (the plugin's runner
    tests); `pack lint` refuses it in a pack meant for others, and any other http address is refused outright."""
    monkeypatch.setattr(E, "source_url", REAL_SOURCE_URL)
    pack = remote.make(remote.planes(base_url=remote.base))
    assert pack.problems == [] and E.is_loopback_url(remote.base)
    assert packs.fetch_example(pack, "planes", remote.out)["status"] == "fetched"
    assert any("is not https" in p for p in packs.lint_pack(pack.folder)["problems"])
    assert not E.is_loopback_url("http://data.example.invalid/") and not E.is_loopback_url("https://127.0.0.1/")


def test_check_expect(tmp_path):
    spec = {"expect": {"n_objects": [2, 6], "spots_per_object": [1.0, 3.0]}}
    (tmp_path / "x__image_summary.csv").write_text("n_nuclei,mean_spots,qc_pass\n4,2.5,True\n", encoding="utf-8")
    assert E.check_expect(spec, tmp_path)["ok"]
    (tmp_path / "x__image_summary.csv").write_text("n_nuclei,mean_spots,qc_pass\n9,2.5,True\n", encoding="utf-8")
    res = E.check_expect(spec, tmp_path)
    assert not res["ok"] and [c["key"] for c in res["checks"] if not c["ok"]] == ["n_objects"]
    assert E.check_expect({}, tmp_path) == {"ok": True, "checks": []}


@pytest.mark.parametrize("change, needle", [
    (lambda s: s.update(base_url="http://data.example.invalid/"), "must be an https:// address"),
    (lambda s: s["sources"][0].pop("sha256"), "sha256 must be 64"),
    (lambda s: s["sources"][0].update(size=0), "size must be the file's size"),
    (lambda s: s["sources"][0].update(path="../escape.tif"), "must be a relative path"),
    (lambda s: s["sources"][0].update(crc32="abc"), "crc32 must be 8 lowercase hex digits"),
    (lambda s: s["sources"][0].update(zip_size=10), "crc32 and zip_size belong to a zip member"),
    (lambda s: s["assemble"].update(axes="ZQYX"), "assemble.axes"),
    (lambda s: s["assemble"].update(channels=["a", "a"]), "unique channel names"),
    (lambda s: s["assemble"].update(spacing_zyx_um=[0.5, 0.2]), "three positive numbers"),
    (lambda s: s["assemble"].update(crop={"y": [5, 2]}), "assemble.crop y"),
    (lambda s: s.update(pixels_sha256="nope"), "pixels_sha256 must be 64"),
    (lambda s: s.pop("licence"), "licence is missing"),
    (lambda s: s.pop("attribution"), "attribution is missing"),
    (lambda s: s.update(expect={"n_objects": [5, 1]}), "expect.n_objects must be [low, high]"),
    (lambda s: s.update(expect={"cells": [1, 5]}), "expect.cells: unknown"),
    (lambda s: s.update(recipe="nope"), "recipe 'nope' is not a recipe of the pack"),
    (lambda s: s.update(image="example/x.ome.tif"), "either image"),
    (lambda s: s.update(id="Bad-Id"), "id does not match"),
])
def test_malformed_example_specs_make_the_pack_unusable(remote, change, needle):
    spec = remote.planes()
    change(spec)
    data = yaml.safe_load((TEST_PACK / "pack.yaml").read_text(encoding="ascii"))
    probs = packs.validate_manifest({**data, "examples": [*data["examples"], spec]})
    assert any(needle in p for p in probs), probs
    assert packs.validate_manifest({**data, "examples": [*data["examples"], remote.planes()]}) == []
    dup = {**data, "examples": [*data["examples"], remote.planes(), remote.planes()]}
    assert any("declared twice" in p for p in packs.validate_manifest(dup))


def test_shipped_example_specs_are_valid_and_sources_are_named_by_file():
    for pid in ("cultured_cells", "budding_yeast"):
        pack = packs.load_pack(ROOT / "packs" / pid)
        remote = [e for e in packs.list_examples(pack) if e["kind"] == "remote"]
        assert len(remote) == 1 and pack.problems == []
        for src in remote[0]["sources"]:
            assert E.source_name(src) == (src.get("member") or src["path"]).split("/")[-1]
            assert E.source_url(remote[0], src).startswith("https://")
    assert E.source_name({"path": "foci_detection_validation_subset.zip/content",
                          "member": "valid/data_258.tif"}) == "data_258.tif"
