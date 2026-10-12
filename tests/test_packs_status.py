"""pack status, `pack check` and the pack check of `run --pack` (design 3.9, P1 acceptance 6).

pack status reports a missing module (find_spec), an undeclared backend, a missing readers extra, a missing
model, an engine out of range, a recipe hash mismatch, a faulty tunable, an engine file that differs from its
recorded digest and a pixel size outside expects (a warning), each with its fix; the image header is read
without bioio (an .nd2 through the nd2 package, a TIFF through tifffile, anything else only by the deep check);
the deep check reports the device and warns when an accelerator is recommended and absent; `pack check` exits 2
when unusable or unknown; `run --pack` refuses an unusable pack unless --no-pack-check.
"""
from __future__ import annotations

import json

import numpy as np
import pytest
import tifffile
from packs_helpers import (  # noqa: F401  (pack_env is a fixture)
    CELEGANS_PACK,
    TEST_PACK,
    copy_pack,
    pack_env,
    patch_nd2_reader,
    small_stack,
    use_fixture_packs,
)

from germquant import cli, cli_packs, packs
from germquant.io import ome_header as OH
from germquant.packs import status as S

NAME = "20261011_field_001.nd2"


def _env_pack(monkeypatch, tmp_path, edit, rebuild=True, version="0.1.0"):
    root = tmp_path / "env_packs"
    copy_pack(TEST_PACK, root / "test_pack" / version, edit=edit, rebuild=rebuild)
    monkeypatch.setenv("GERMQUANT_PACKS", str(root))
    packs.refresh()
    return root / "test_pack" / version


def _whats(st) -> str:
    return " | ".join(b["what"] + " -> " + b["fix"] for b in st["blockers"])


def test_basic_recipe_is_usable_anywhere(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    st = packs.pack_status("test_pack")
    assert st["usable"] and st["recipe"] == "basic" and st["blockers"] == [] and st["warnings"] == []


def test_missing_module_blocks_with_the_extra_to_install(pack_env, monkeypatch, tmp_path):  # noqa: F811
    _env_pack(monkeypatch, tmp_path, lambda d: d["recipes"][0]["needs"].update(modules=["skan", "no_such_module_xyz"]))
    st = packs.pack_status("test_pack")
    assert not st["usable"] and "no_such_module_xyz" in _whats(st)
    assert "skan" not in _whats(st) or "'sc'" in _whats(st)


def test_missing_model_blocks_with_fetch_models(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    st = packs.pack_status("test_pack", recipe="with_model")
    assert not st["usable"] and st["models"][0]["state"] == "missing"
    assert any(b["fix"] == "germquant pack fetch-models test_pack --recipe with_model" for b in st["blockers"])


def test_engine_out_of_range_blocks(pack_env, monkeypatch, tmp_path):  # noqa: F811
    _env_pack(monkeypatch, tmp_path, lambda d: d.update(requires={"germquant": ">=9.0,<10"}))
    st = packs.pack_status("test_pack")
    assert not st["usable"] and "outside the pack's range >=9.0,<10" in _whats(st)


def test_recipe_hash_mismatch_blocks(pack_env, monkeypatch, tmp_path):  # noqa: F811
    _env_pack(monkeypatch, tmp_path, lambda d: d["recipes"][0].update(expected_config_hash="0123456789ab"))
    st = packs.pack_status("test_pack")
    assert not st["usable"] and "the pack expects 0123456789ab" in _whats(st)


def test_changed_file_blocks_through_the_sums(pack_env, monkeypatch, tmp_path):  # noqa: F811
    folder = _env_pack(monkeypatch, tmp_path, None)
    (folder / "README.md").write_bytes(b"edited after the sums were built\n")
    st = packs.pack_status("test_pack")
    assert not st["usable"] and "SHA256SUMS" in _whats(st)


def test_pixel_size_outside_expects_is_a_warning(pack_env, monkeypatch, tmp_path):  # noqa: F811
    import germquant.io as gio

    use_fixture_packs(monkeypatch)
    meta = {"spacing": (0.5, 0.65, 0.65), "spacing_ok": True, "channel_names": ["DAPI", "foci"]}
    monkeypatch.setattr(gio, "read_image_metadata", lambda path, **k: meta)
    image = tmp_path / NAME
    image.write_bytes(b"x")
    st = packs.pack_status("test_pack", image=image)
    assert st["usable"] and len(st["warnings"]) == 1
    assert "pixel size 0.65 um lies outside" in st["warnings"][0]["what"]
    meta["spacing"] = (0.5, 0.2, 0.2)
    assert packs.pack_status("test_pack", image=image)["warnings"] == []


def test_pack_check_exit_codes(pack_env, monkeypatch, capsys):  # noqa: F811
    use_fixture_packs(monkeypatch)
    assert cli_packs.main(["check", "test_pack"]) == 0
    assert cli_packs.main(["check", "test_pack", "--recipe", "with_model"]) == 2
    out = capsys.readouterr().out
    assert "NOT usable" in out and "fix: germquant pack fetch-models test_pack" in out
    assert cli.main(["pack", "check", "test_pack@0.1.0", "--recipe", "with_model", "--json"]) == 2
    capsys.readouterr()
    assert cli_packs.main(["check", "nope"]) == 2                        # an unknown pack is unusable too
    assert "no pack 'nope'" in capsys.readouterr().err
    assert cli.main(["pack", "check", "nope", "--json"]) == 2
    st = json.loads(capsys.readouterr().out)
    assert st["usable"] is False and "no pack 'nope'" in st["blockers"][0]["what"]


def test_run_pack_refuses_an_unusable_pack_unless_no_pack_check(pack_env, monkeypatch, tmp_path, capsys):  # noqa: F811
    _env_pack(monkeypatch, tmp_path, lambda d: d.update(requires={"germquant": ">=9.0,<10"}))
    patch_nd2_reader(monkeypatch, small_stack())
    with pytest.raises(SystemExit) as e:
        cli.main(["run", NAME, "--pack", "test_pack", "--out", str(tmp_path / "refused")])
    assert e.value.code == 2 and not (tmp_path / "refused").exists()
    err = capsys.readouterr().err
    assert "REFUSED" in err and "outside the pack's range" in err and "--no-pack-check" in err
    assert cli.main(["run", NAME, "--pack", "test_pack", "--no-pack-check", "--out", str(tmp_path / "ran")]) == 0
    assert (tmp_path / "ran" / "run_manifest.json").is_file()


def test_recipe_needs_pack(pack_env, capsys):  # noqa: F811
    with pytest.raises(SystemExit) as e:
        cli.main(["run", NAME, "--config", "config/config.yaml", "--recipe", "basic", "--out", "x"])
    assert e.value.code == 2 and "--recipe needs --pack" in capsys.readouterr().err


# ---- backends, the readers extra, tunables and engine files ----------------------------------------------------

def test_undeclared_backends_block_with_the_entry_point_fix(pack_env, monkeypatch, tmp_path):  # noqa: F811
    _env_pack(monkeypatch, tmp_path, lambda d: d["recipes"][0]["needs"].update(
        segment_backends=["yeaz_not_installed"], spot_backends=["spotiflow_not_installed"]))
    st = packs.pack_status("test_pack")
    text = _whats(st)
    assert not st["usable"]
    assert "segmentation backend 'yeaz_not_installed' is not installed" in text
    assert "an entry point in germquant.segment_backends" in text
    assert "spot backend 'spotiflow_not_installed' is not installed" in text
    assert "an entry point in germquant.spot_backends" in text


def test_a_non_nd2_image_needs_the_readers_extra(pack_env, monkeypatch, tmp_path):  # noqa: F811
    use_fixture_packs(monkeypatch)
    monkeypatch.setattr(S, "_module_present", lambda name: name != "bioio")
    st = packs.pack_status("test_pack", image=tmp_path / "field.ome.tif")
    assert not st["usable"] and "is not an .nd2 and the readers extra is not installed" in _whats(st)
    assert "install the 'readers' extra" in _whats(st) and st["warnings"] == []
    assert packs.pack_status("test_pack", image=tmp_path / NAME)["usable"]       # an .nd2 needs no extra


def test_faulty_tunables_block_the_recipe(pack_env, monkeypatch, tmp_path):  # noqa: F811
    _env_pack(monkeypatch, tmp_path, lambda d: d["recipes"][0].update(
        tunables=[{"key": "spots.no_such_key", "label": "typo"}, {"key": "spots.log", "label": "a section"}]))
    st = packs.pack_status("test_pack")
    text = _whats(st)
    assert not st["usable"] and "use a corrected pack version" in text
    assert "recipe basic: tunable spots.no_such_key is not set by the recipe" in text
    assert "recipe basic: tunable spots.log is not a scalar" in text


def test_engine_file_digest_mismatch_is_reported_everywhere(pack_env, monkeypatch, tmp_path):  # noqa: F811
    from germquant.packs.install import PackInstallError

    def edit(d):
        d["version"] = "0.1.1"                                        # not the built-in id and version
        d["engine_files"][sorted(d["engine_files"])[0]] = "0" * 64
    root = tmp_path / "env_packs"
    folder = copy_pack(CELEGANS_PACK, root / "celegans_germline" / "0.1.1", edit=edit)
    monkeypatch.setenv("GERMQUANT_PACKS", str(root))
    packs.refresh()
    needle = "differs from the digest the pack records"
    assert any(needle in p for p in packs.verify_pack(folder))
    st = packs.pack_status("celegans_germline@0.1.1")
    assert not st["usable"] and needle in _whats(st)
    assert any(needle in p for p in packs.lint_pack(folder)["problems"])
    with pytest.raises(PackInstallError, match=needle):
        packs.vendor("celegans_germline", "0.1.1", tmp_path / "project" / "packs")
    assert not (tmp_path / "project" / "packs" / "celegans_germline").exists()
    assert packs.engine_files_problems(packs.find_pack("celegans_germline@0.1.0")) == []   # the built-in is fine


# ---- the image header without bioio, and the deep check ---------------------------------------------------------

def _ome_tiff(path, spacing, names=("DAPI", "foci")):
    data = np.zeros((len(names), 3, 16, 16), np.uint16)
    tifffile.imwrite(path, data, ome=True, metadata={
        "axes": "CZYX", "PhysicalSizeZ": spacing[0], "PhysicalSizeY": spacing[1], "PhysicalSizeX": spacing[2],
        "Channel": {"Name": list(names)}})
    return path


def _imagej_tiff(path, spacing, unit="micron"):
    data = np.zeros((3, 2, 16, 16), np.uint16)                                      # ZCYX
    tifffile.imwrite(path, data, imagej=True, resolution=(1 / spacing[2], 1 / spacing[1]),
                     metadata={"axes": "ZCYX", "spacing": spacing[0], "unit": unit})
    return path


def _plain_tiff(path):
    tifffile.imwrite(path, np.zeros((16, 16), np.uint16))
    return path


def test_tiff_headers_are_checked_without_bioio(pack_env, monkeypatch, tmp_path):  # noqa: F811
    import germquant.io as gio

    use_fixture_packs(monkeypatch)
    monkeypatch.setattr(S, "_module_present", lambda name: True)    # as where the readers extra is installed

    def through_bioio(*a, **k):
        raise AssertionError("the shallow check read a TIFF header through bioio")

    monkeypatch.setattr(gio, "read_image_metadata", through_bioio)
    for ok in (_ome_tiff(tmp_path / "ok.ome.tif", (0.5, 0.2, 0.2)), TEST_PACK / "example" / "test_pack_synthetic.ome.tif"):
        assert packs.pack_status("test_pack", image=ok)["warnings"] == [], ok
    coarse = _ome_tiff(tmp_path / "coarse.tif", (0.5, 0.65, 0.65))                # OME-XML under a plain .tif name
    st = packs.pack_status("test_pack", image=coarse)
    assert st["usable"] and len(st["warnings"]) == 1 and "pixel size 0.65 um lies outside" in st["warnings"][0]["what"]
    st = packs.pack_status("test_pack", image=_imagej_tiff(tmp_path / "ij.tif", (2.5, 0.2, 0.2)))
    assert len(st["warnings"]) == 1 and "z step 2.5 um lies outside" in st["warnings"][0]["what"]
    st = packs.pack_status("test_pack", image=_imagej_tiff(tmp_path / "ij_nm.tif", (500.0, 650.0, 650.0), unit="nm"))
    assert len(st["warnings"]) == 1 and "pixel size 0.65 um lies outside" in st["warnings"][0]["what"]
    assert packs.pack_status("test_pack", image=_plain_tiff(tmp_path / "plain.tif"))["warnings"] == []   # no size
    st = packs.pack_status("test_pack", image=tmp_path / "missing.tif")
    assert st["usable"] and "could not read the header" in st["warnings"][0]["what"]


@pytest.mark.parametrize("kind", ["ome", "ome_in_plain_tif", "imagej", "imagej_nm", "plain"])
def test_light_tiff_header_matches_the_bioio_reader(tmp_path, kind):
    pytest.importorskip("bioio")
    pytest.importorskip("bioio_ome_tiff")
    pytest.importorskip("bioio_tifffile")
    from germquant.io import read_image_metadata

    path = {"ome": lambda: _ome_tiff(tmp_path / "a.ome.tif", (0.3, 0.11, 0.11)),
            "ome_in_plain_tif": lambda: _ome_tiff(tmp_path / "b.tif", (0.45, 0.07, 0.07)),
            "imagej": lambda: _imagej_tiff(tmp_path / "c.tif", (0.4, 0.2, 0.2)),
            "imagej_nm": lambda: _imagej_tiff(tmp_path / "d.tif", (300.0, 90.0, 90.0), unit="nm"),
            "plain": lambda: _plain_tiff(tmp_path / "e.tif")}[kind]()
    light, full = OH.read_tiff_header(path), read_image_metadata(path)
    assert light["spacing_ok"] == full["spacing_ok"]
    assert light["spacing"] == pytest.approx(full["spacing"], rel=1e-9)
    if kind.startswith("ome"):
        assert light["channel_names"] == full["channel_names"] == ["DAPI", "foci"]


def test_other_formats_are_left_to_the_deep_check(pack_env, monkeypatch, tmp_path):  # noqa: F811
    import germquant.io as gio

    use_fixture_packs(monkeypatch)
    monkeypatch.setattr(S, "_module_present", lambda name: True)
    seen = []

    def header(path, **k):
        seen.append(path)
        return {"spacing": (0.5, 0.65, 0.65), "spacing_ok": True}

    monkeypatch.setattr(gio, "read_image_metadata", header)
    image = tmp_path / "field.czi"
    image.write_bytes(b"x")
    st = packs.pack_status("test_pack", image=image)
    assert st["usable"] and seen == [] and len(st["warnings"]) == 1          # bioio would load dask (and CuPy)
    assert "were not checked" in st["warnings"][0]["what"] and "--deep" in st["warnings"][0]["fix"]
    deep = packs.pack_status("test_pack", image=image, deep=True)
    assert seen == [image] and len(deep["warnings"]) == 1                    # read after the pipeline (torch)
    assert "pixel size 0.65 um lies outside" in deep["warnings"][0]["what"]
    assert deep["device"] in ("cpu", "cuda", "mps")
    # a pack that expects nothing reads no header at all

    def no_expects(d):
        d.pop("expects")
        d["version"] = "0.3.0"

    _env_pack(monkeypatch, tmp_path, no_expects, version="0.3.0")
    assert packs.pack_status("test_pack@0.3.0", image=image, deep=True)["warnings"] == [] and seen == [image]


def test_deep_check_reports_the_device_and_warns_without_an_accelerator(pack_env, monkeypatch, tmp_path):  # noqa: F811
    from germquant import device

    _env_pack(monkeypatch, tmp_path, lambda d: d["recipes"][0]["needs"].update(accelerator="recommended"))
    monkeypatch.setattr(device, "available", lambda: {"cuda": False, "mps": False})
    shallow = packs.pack_status("test_pack")
    assert shallow["warnings"] == [] and shallow["device"] is None               # the shallow check never looks
    st = packs.pack_status("test_pack", deep=True)
    assert st["usable"] and st["device"] == "cpu"
    assert any("recommends an accelerator" in w["what"] for w in st["warnings"])
    monkeypatch.setattr(device, "available", lambda: {"cuda": True, "mps": False})
    st = packs.pack_status("test_pack", deep=True)
    assert st["device"] == "cuda" and st["warnings"] == []
