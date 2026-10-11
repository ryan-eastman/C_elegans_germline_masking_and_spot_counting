"""The two extends fallbacks, the recorded anchors and the portable stage payload (design 3.4, P1 acceptance 3).

* The test pack copied under two roots loads with one config_hash.
* ``extends: generic.yaml`` from outside the checkout resolves to the checkout file (the bare-name fallback);
  ``extends: pack:test_pack@0.1.0/basic`` resolves; an unknown pack or recipe raises FileNotFoundError that
  lists the candidates tried and the packs available.
* Every config and profile of the checkout resolves its extends to the same file as today's rule (rule R4).
* An edited recipe extending pack:celegans_germline@0.1.0/default gives the same stage hashes when the engine
  config root moves (the built-in root monkeypatched to one temporary copy of config/ and the pack, then to
  another; a fixed model digest); the same for an edited recipe on the test pack copied under two roots.
* An out-of-repo overlay with an absolute extends path keeps today's stage hashes.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest
import yaml
from packs_helpers import CELEGANS_PACK, ROOT, TEST_PACK, copy_pack, pack_env  # noqa: F401

from germquant import config as C
from germquant import packs, stages
from germquant.config import load_config

FIXED = {tool: "PINNED" for tools in stages._STAGE_TOOLS.values() for tool in tools}
ROLES = {"dna": 0, "central_element": 1, "foci": 2, "granule": None, "lamin": None}


def _hashes(cfg) -> dict:
    return stages.stage_hashes(cfg, ROLES, FIXED, model_sha="fixed-model-digest")


def _legacy_resolve(path: Path, base_dir: Path, ext: str) -> Path:
    """Today's extends rule, copied from the base commit (bc62476), to compare with."""
    ext_p = Path(ext)
    if ext_p.is_absolute():
        cand = [ext_p]
    elif len(ext_p.parts) > 1:
        cand = [path.parent / ext_p, base_dir / "config" / ext_p, base_dir / ext_p]
    else:
        cand = [base_dir / "config" / ext_p, base_dir / ext_p, path.parent / ext_p]
    for c in cand:
        if c.is_file() and c.resolve() != path.resolve():
            return c
    raise FileNotFoundError(ext)


def test_every_checkout_config_resolves_to_the_same_file_as_before():
    files = sorted((ROOT / "config").glob("*.yaml")) + sorted((ROOT / "config" / "profiles").glob("*.yaml"))
    checked = 0
    for f in files:
        ext = (yaml.safe_load(f.read_text(encoding="cp1252")) or {}).get("extends")
        if not ext:
            continue
        got, new_form = C._resolve_extends_ex(f, C._base_dir_for(f), str(ext))
        assert (got, new_form) == (_legacy_resolve(f, C._base_dir_for(f), str(ext)), False), f
        assert load_config(f).anchors == {}, f                 # nothing anchored, today's payload
        checked += 1
    assert checked >= 7


def test_pack_copied_under_two_roots_loads_with_one_hash(tmp_path):
    a = copy_pack(TEST_PACK, tmp_path / "rootA" / "test_pack", rebuild=False)
    b = copy_pack(TEST_PACK, tmp_path / "somewhere else" / "rootB" / "test_pack", rebuild=False)
    for rid in ("basic", "with_model"):
        ha = load_config(a / "config" / "profiles" / f"{rid}.yaml")
        hb = load_config(b / "config" / "profiles" / f"{rid}.yaml")
        assert ha.hash == hb.hash and _hashes(ha) == _hashes(hb)
        assert ha.base_dir == a and ha.anchors == {}


def test_bare_name_resolves_to_the_checkout_from_anywhere(pack_env, tmp_path):  # noqa: F811
    study = tmp_path / "study"
    study.mkdir()
    prof = study / "x.yaml"
    prof.write_text("extends: generic.yaml\nspots:\n  effect_size_min: 4.0\n", encoding="ascii")
    got, new_form = C._resolve_extends_ex(prof, C._base_dir_for(prof), "generic.yaml")
    assert got == ROOT.resolve() / "config" / "profiles" / "generic.yaml" and new_form is True
    cfg = load_config(prof)
    gen = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    assert cfg.get("spots.effect_size_min") == 4.0 and cfg.get("germline.enabled") is False
    assert cfg.anchors["portable"] is True and "io.channel_map" in cfg.anchors["paths"]
    assert cfg.anchors["paths"]["io.channel_map"]["relative"] == gen.get("io.channel_map")
    # the read stage hashes the relative channel-map path, so the chain hashes alike on every machine
    sec = cfg.get("io")
    assert Path(sec["channel_map"]).is_absolute()


def test_pack_form_resolves_and_unknown_names_list_candidates(pack_env, monkeypatch, tmp_path):  # noqa: F811
    monkeypatch.setenv("GERMQUANT_PACKS", str(TEST_PACK.parent))
    packs.refresh()
    prof = tmp_path / "edit.yaml"
    prof.write_text("extends: pack:test_pack@0.1.0/basic\nspots:\n  log:\n    min_snr: 3.0\n", encoding="ascii")
    cfg = load_config(prof)
    assert cfg.get("spots.log.min_snr") == 3.0 and cfg.get("segmentation.nuclei.method") == "classical"
    base = load_config(TEST_PACK / "config" / "profiles" / "basic.yaml")
    assert cfg.hash != base.hash and cfg.anchors["portable"] is True
    for bad, needle in (("pack:nope@0.1.0/basic", r"Packs available: .*test_pack@0\.1\.0"),
                        ("pack:test_pack@0.1.0/nope", "its recipes: basic, with_model"),
                        ("nothing_like_this.yaml", r"tried \[.*nothing_like_this.yaml.*packs available")):
        p = tmp_path / "bad.yaml"
        p.write_text(f"extends: {bad}\n", encoding="ascii")
        with pytest.raises(FileNotFoundError, match=needle):
            load_config(p)


def _engine_root_copy(dest: Path) -> Path:
    shutil.copytree(ROOT / "config", dest / "config")
    shutil.copytree(CELEGANS_PACK, dest / "packs" / "celegans_germline")
    (dest / "src" / "germquant").mkdir(parents=True)
    return dest


def test_edited_celegans_recipe_hashes_alike_when_the_engine_root_moves(pack_env, monkeypatch, tmp_path):  # noqa: F811
    study = tmp_path / "study"
    study.mkdir()
    prof = study / "edited.yaml"
    prof.write_text("extends: pack:celegans_germline@0.1.0/default\nspots:\n  effect_size_min: 3.5\n", encoding="ascii")
    results = []
    for name in ("engine_a", "deeper/engine b"):
        root = _engine_root_copy(tmp_path / name)
        monkeypatch.setattr(C, "builtin_roots", lambda root=root: [root.resolve()])
        packs.refresh()
        cfg = load_config(prof)
        assert cfg.get("spots.effect_size_min") == 3.5
        assert str(root.resolve()) in cfg.get("io.channel_map")       # the absolute path differs per root ...
        results.append((cfg.hash, _hashes(cfg), cfg.get("io.channel_map"),
                        cfg.get("segmentation.nuclei.cellpose_model")))
    (ha, sa, cma, ma), (hb, sb, cmb, mb) = results
    assert cma != cmb and ma != mb
    assert ha == hb                                                    # ... the hashes do not
    assert sa == sb, {k: (sa[k], sb[k]) for k in sa if sa[k] != sb[k]}


def test_edited_test_pack_recipe_hashes_alike_under_two_roots(pack_env, monkeypatch, tmp_path):  # noqa: F811
    prof = tmp_path / "edited.yaml"
    prof.write_text("extends: pack:test_pack@0.1.0/basic\nsegmentation:\n  nuclei:\n    min_volume_um3: 3.0\n",
                    encoding="ascii")
    seen = []
    for name in ("pa", "pb with space"):
        root = tmp_path / name
        copy_pack(TEST_PACK, root / "test_pack" / "0.1.0")
        monkeypatch.setenv("GERMQUANT_PACKS", str(root))
        packs.refresh()
        cfg = load_config(prof)
        seen.append((cfg.hash, _hashes(cfg), cfg.get("io.channel_map")))
    assert seen[0][2] != seen[1][2]
    assert seen[0][0] == seen[1][0] and seen[0][1] == seen[1][1]


def test_absolute_extends_overlay_keeps_todays_stage_hashes(tmp_path):
    prof = tmp_path / "study.yaml"
    prof.write_text(f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}\ncoloc:\n  enabled: false\n",
                    encoding="ascii")
    cfg = load_config(prof)
    assert cfg.anchors["portable"] is False and cfg.anchors["paths"]          # anchored, recorded, not portable
    today = _hashes(cfg)
    cfg.anchors = {}                                                         # what the base commit had
    assert _hashes(cfg) == today
    assert Path(cfg.get("io.channel_map")).is_absolute()


def test_portable_payload_respects_overrides(pack_env, tmp_path):  # noqa: F811
    prof = tmp_path / "x.yaml"
    prof.write_text("extends: generic.yaml\n", encoding="ascii")
    cfg = load_config(prof)
    before = _hashes(cfg)
    cfg.set("io.channel_map", str(ROOT / "config" / "channel_maps" / "n2_dapi_syp3_rad51.yaml"))
    assert _hashes(cfg)["read"] != before["read"]                 # an overridden path keeps its own value
