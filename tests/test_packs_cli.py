"""The command surface of packs (design 3.13 and 5.1): `germquant pack ...`, the light `python -m germquant.packs`
entry, and declared_backend_names on both registries, which lists entry points without loading them."""
from __future__ import annotations

import json

import pytest
from packs_helpers import TEST_PACK, TINYNET, TINYNET_BYTES, copy_pack, pack_env, use_fixture_packs  # noqa: F401

from germquant import cli, cli_packs, packs
from germquant.segment import backends as SB
from germquant.spots import backends as SPB


class _EP:
    def __init__(self, name):
        self.name, self.value = name, f"pkg:{name}"

    def load(self):                                    # pragma: no cover - must never be called
        raise AssertionError("declared_backend_names loaded an entry point")


def test_declared_backend_names_do_not_load_entry_points(monkeypatch):
    monkeypatch.setattr(SB, "_iter_entry_points", lambda: [_EP("yeaz"), _EP("cellpose")])
    monkeypatch.setattr(SPB, "_iter_entry_points", lambda: [_EP("spotiflow"), _EP("bad name!")])
    seg = SB.declared_backend_names()
    assert {"cellpose", "classical", "nuclei3d", "yeaz"} <= set(seg) and seg == sorted(seg)
    spots = SPB.declared_backend_names()
    assert {"spotmax", "log", "net3d", "spotiflow"} <= set(spots) and "bad name!" not in spots


def test_pack_commands_end_to_end(pack_env, monkeypatch, tmp_path, capsys):  # noqa: F811
    use_fixture_packs(monkeypatch)
    assert cli.main(["pack", "list"]) == 0
    out = capsys.readouterr().out
    assert "test_pack 0.1.0" in out and "celegans_germline 0.1.0" in out and "recipes: basic*, with_model" in out
    assert cli.main(["pack", "list", "--json", "--all"]) == 0
    listing = json.loads(capsys.readouterr().out)
    assert listing["builtin_root_found"] and listing["version"]["germquant"]
    assert cli.main(["pack", "show", "test_pack", "--recipe", "with_model"]) == 0
    assert "model tinynet@f5055f8a" in capsys.readouterr().out
    assert cli.main(["pack", "verify", "test_pack"]) == 0
    assert cli.main(["pack", "lint", str(TEST_PACK)]) == 0
    # fetch-models with --from-file, then the recipe is usable
    f = tmp_path / "tinynet.bin"
    f.write_bytes(TINYNET_BYTES)
    assert cli.main(["pack", "fetch-models", "test_pack", "--recipe", "with_model", "--no-download"]) == 1
    assert cli.main(["pack", "fetch-models", "test_pack", "--recipe", "with_model", "--from-file", f"{TINYNET}={f}"]) == 0
    capsys.readouterr()
    assert cli.main(["pack", "check", "test_pack", "--recipe", "with_model"]) == 0
    # build-sums, install into the first $GERMQUANT_PACKS folder, refused second install with other content
    work = copy_pack(TEST_PACK, tmp_path / "work" / "test_pack", edit=lambda d: d.update(version="0.1.5"))
    assert cli.main(["pack", "build-sums", str(work)]) == 0
    target = tmp_path / "installed"
    monkeypatch.setenv("GERMQUANT_PACKS", str(target))
    packs.refresh()
    assert cli.main(["pack", "install", str(work)]) == 0
    assert (target / "test_pack" / "0.1.5" / "pack.yaml").is_file()
    (work / "README.md").write_bytes(b"x\n")
    assert cli.main(["pack", "build-sums", str(work)]) == 0
    assert cli.main(["pack", "install", str(work)]) == 1
    assert "immutable" in capsys.readouterr().err


def test_light_entry_has_no_golden_and_prints_versions(capsys):
    assert cli_packs.main(["--version"]) == 0
    assert capsys.readouterr().out.startswith("germquant ")
    with pytest.raises(SystemExit):
        cli_packs.main(["golden", "test_pack"])
