"""pack status, `pack check` and the pack check of `run --pack` (design 3.9, P1 acceptance 6).

pack status reports a missing module (find_spec), a missing model, an engine out of range, a recipe hash
mismatch and a pixel size outside expects (a warning), each with its fix; `pack check` exits 2 when unusable;
`run --pack` refuses an unusable pack unless --no-pack-check.
"""
from __future__ import annotations

import pytest
from packs_helpers import TEST_PACK, copy_pack, pack_env, patch_nd2_reader, small_stack, use_fixture_packs  # noqa: F401

from germquant import cli, cli_packs, packs

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
    assert cli_packs.main(["check", "nope"]) == 1                        # an unknown pack: a clean error


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
