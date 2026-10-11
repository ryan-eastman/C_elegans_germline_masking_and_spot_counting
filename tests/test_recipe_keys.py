"""germquant.recipe_keys: the strict checks of the opt-in recipe keys of tranche 5 package P2a
(docs/DESIGN_system_packs.md sections 3.10, 3.12, 4.2, 4.4 and 4.7; docs/PLANAR.md, "Recipe keys for packs").

* No shipped config or profile sets any of the keys, so a default run checks nothing.
* Valid values (true, false, null, numbers, numpy scalars) pass; every malformed value is refused with its key
  named, a quoted "false" or "4" included, and every problem is listed at once, spots.instances entries too.
* process_image refuses a malformed key before any stage runs and before any file is written (a stale done
  marker is still cleared). This holds under method auto, where the Cellpose backend's own error would be
  swallowed by the watershed fallback, in a 3D run and in a projected run. A valid object diameter still
  reaches Cellpose under auto.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel synthetic gonad)
from test_pack_knobs import fake_cellpose  # noqa: F401  (fixture: Cellpose mocked)
from test_projected_run import IMAGE, patch_reader, recipe, yeast_like_stack

from germquant import pipeline, recipe_keys
from germquant.config import load_config
from germquant.spots.detect import SPOTMAX_OPTIONAL_KEYS

ROOT = Path(__file__).resolve().parents[1]
ND2 = "20251105_n2_nohs_HERM_001.nd2"
ND2_ID = "20251105_n2_nohs_HERM_001"
SHIPPED = sorted((ROOT / "config").glob("*.yaml")) + sorted((ROOT / "config" / "profiles").glob("*.yaml"))


class _Cfg:
    """A stand-in config over a nested mapping, with Config.get's dotted lookup."""

    def __init__(self, data: dict):
        self._data = data

    def get(self, dotted, default=None):
        node = self._data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node


def _with(dotted: str, value) -> _Cfg:
    data: dict = {}
    node = data
    *parents, last = dotted.split(".")
    for part in parents:
        node = node.setdefault(part, {})
    node[last] = value
    return _Cfg(data)


@pytest.mark.parametrize("path", SHIPPED, ids=lambda p: p.name)
def test_no_shipped_config_sets_a_checked_key(path):
    cfg = load_config(path)
    assert [k for k in recipe_keys.CHECKED_KEYS if cfg.get(k) is not None] == []
    for inst in cfg.get("spots.instances", []) or []:
        assert not set(inst) & set(SPOTMAX_OPTIONAL_KEYS)
    assert recipe_keys.problems(cfg) == []


def test_the_spot_keys_are_the_optional_spotmax_keys():
    assert set(recipe_keys.SPOT_LENGTH_KEYS + recipe_keys.SPOT_FLAG_KEYS) == set(SPOTMAX_OPTIONAL_KEYS)


GOOD = [
    ("segmentation.nuclei.cellpose", {}), ("segmentation.nuclei.cellpose", {"object_diameter_um": None}),
    ("segmentation.nuclei.cellpose.object_diameter_um", 12), ("segmentation.nuclei.cellpose.object_diameter_um", 0.5),
    ("segmentation.nuclei.cellpose.object_diameter_um", np.float64(4.0)),    # a float subclass, JSON-safe
    ("segmentation.nuclei.require_channel", True), ("segmentation.nuclei.require_channel", False),
    ("qc.min_nuclei", 0), ("qc.min_nuclei", 1), ("qc.min_nuclei", 10),
    ("spots.spot_radius_z_um", 1.0), ("spots.spot_radius_z_um", 2), ("spots.remove_hot_pixels", True),
    ("spots.sharpen", False), ("spots.sharpen", True),
    *((key, None) for key in recipe_keys.CHECKED_KEYS),
]


@pytest.mark.parametrize("key,value", GOOD)
def test_valid_values_pass(key, value):
    cfg = _with(key, value)
    assert recipe_keys.problems(cfg) == []
    recipe_keys.check(cfg)


BAD = [
    ("segmentation.nuclei.cellpose", 4.0), ("segmentation.nuclei.cellpose", "big"),
    ("segmentation.nuclei.cellpose.object_diameter_um", 0), ("segmentation.nuclei.cellpose.object_diameter_um", -5.0),
    ("segmentation.nuclei.cellpose.object_diameter_um", "ten"), ("segmentation.nuclei.cellpose.object_diameter_um", "4"),
    ("segmentation.nuclei.cellpose.object_diameter_um", True),
    ("segmentation.nuclei.cellpose.object_diameter_um", float("inf")),
    ("segmentation.nuclei.cellpose.object_diameter_um", float("nan")),
    ("segmentation.nuclei.require_channel", "false"), ("segmentation.nuclei.require_channel", "true"),
    ("segmentation.nuclei.require_channel", 1), ("segmentation.nuclei.require_channel", "yes"),
    ("qc.min_nuclei", 2.7), ("qc.min_nuclei", 3.0), ("qc.min_nuclei", -1), ("qc.min_nuclei", True),
    ("qc.min_nuclei", "10"),
    ("spots.spot_radius_z_um", 0), ("spots.spot_radius_z_um", -1.0), ("spots.spot_radius_z_um", "1.0"),
    ("spots.spot_radius_z_um", False),
    ("spots.sharpen", "false"), ("spots.sharpen", 0), ("spots.remove_hot_pixels", "no"),
    ("spots.remove_hot_pixels", 1),
    # numpy scalars other than float64 never come from YAML, and Config.hash (JSON) could not record them
    ("segmentation.nuclei.require_channel", np.bool_(True)), ("qc.min_nuclei", np.int64(3)),
    ("segmentation.nuclei.cellpose.object_diameter_um", np.float32(4.0)),
]


@pytest.mark.parametrize("key,value", BAD)
def test_each_malformed_value_is_refused_with_its_key_named(key, value):
    cfg = _with(key, value)
    (problem,) = recipe_keys.problems(cfg)
    assert problem.startswith(f"{key} must be"), problem
    if isinstance(value, str) and key != "segmentation.nuclei.cellpose":
        assert "a quoted value is text" in problem          # the likely mistake is named
    with pytest.raises(ValueError, match=re.escape(key)):
        recipe_keys.check(cfg)


def test_every_problem_is_listed_at_once_and_instances_are_checked():
    cfg = _Cfg({"segmentation": {"nuclei": {"require_channel": "false", "cellpose": {"object_diameter_um": -1}}},
                "qc": {"min_nuclei": 2.5},
                "spots": {"sharpen": "false",
                          "instances": [{"name": "foci2", "role": "foci2", "sharpen": "no", "spot_radius_z_um": 0},
                                        "not a mapping",                  # pipeline.spot_instances refuses this one
                                        {"name": "ok", "role": "x", "sharpen": False, "spot_radius_z_um": 0.6}]}})
    found = recipe_keys.problems(cfg)
    assert [p.split(" must be")[0] for p in found] == [
        "segmentation.nuclei.cellpose.object_diameter_um", "segmentation.nuclei.require_channel", "qc.min_nuclei",
        "spots.sharpen", "spots.instances[0].spot_radius_z_um", "spots.instances[0].sharpen"]
    with pytest.raises(ValueError, match="recipe refused before any stage ran") as err:
        recipe_keys.check(cfg)
    assert all(p in str(err.value) for p in found)


def test_object_diameter_um():
    for unset in (None, {}, {"cellpose": None}, {"cellpose": {}}, {"cellpose": {"object_diameter_um": None}}):
        assert recipe_keys.object_diameter_um(unset) is None
    assert recipe_keys.object_diameter_um({"cellpose": {"object_diameter_um": 12}}) == 12.0
    with pytest.raises(ValueError, match="segmentation.nuclei.cellpose.object_diameter_um must be"):
        recipe_keys.object_diameter_um({"cellpose": {"object_diameter_um": "12"}})


@pytest.mark.parametrize("key,value", [
    ("segmentation.nuclei.cellpose.object_diameter_um", -2.0),
    ("segmentation.nuclei.cellpose", 4.0),
    ("segmentation.nuclei.require_channel", "false"),
    ("qc.min_nuclei", 2.7),
    ("spots.sharpen", "false"),
])
def test_process_image_refuses_a_malformed_key_before_any_stage(tmp_path, patched_reader, fake_cellpose,  # noqa: F811
                                                                 key, value):
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg.set("segmentation.nuclei.method", "auto")     # the method whose fallback would hide a Cellpose error
    cfg.set("segmentation.nuclei.cellpose_model", "cpsam")
    cfg.set(key, value)
    out = tmp_path / "r"
    out.mkdir()
    (out / f"{ND2_ID}__done.json").write_text("{}", encoding="utf-8")       # a stale marker of an older run
    with pytest.raises(ValueError, match=re.escape(key)):
        pipeline.process_image(ND2, cfg, out, prov=None)
    assert not fake_cellpose.inits and not fake_cellpose.evals
    assert not list(out.iterdir())                    # nothing written, and the stale marker cleared


def test_a_valid_object_diameter_reaches_cellpose_under_auto(tmp_path, patched_reader, fake_cellpose):  # noqa: F811
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg.set("segmentation.nuclei.method", "auto")
    cfg.set("segmentation.nuclei.cellpose_model", "cpsam")
    cfg.set("segmentation.nuclei.cellpose.object_diameter_um", 6.0)
    cfg.set("spots.enabled", False)
    res = pipeline.process_image(ND2, cfg, tmp_path / "r", prov=None)
    assert res["tables"]["image_summary"]["segmentation_method"].iloc[0] == "cellpose"
    (call,) = fake_cellpose.evals
    assert call["kw"] == {"do_3D": True, "z_axis": 0, "channel_axis": None, "anisotropy": 0.4 / 0.2,
                          "diameter": 6.0 / 0.2}


@pytest.mark.parametrize("bad", [{"segmentation.nuclei.cellpose.object_diameter_um": 0},
                                 {"segmentation.nuclei.cellpose.object_diameter_um": -5.0},
                                 {"segmentation.nuclei.cellpose": 4.0}])
def test_a_bad_object_diameter_fails_a_projected_run_under_auto(tmp_path, monkeypatch, fake_cellpose,  # noqa: F811
                                                                 bad):
    stack, _ = yeast_like_stack()
    patch_reader(monkeypatch, stack)
    out = tmp_path / "bad"
    with pytest.raises(ValueError, match="segmentation.nuclei.cellpose"):
        pipeline.process_image(IMAGE, recipe(**{"segmentation.nuclei.method": "auto", "spots.enabled": False, **bad}),
                               out, prov=None)
    assert not fake_cellpose.evals and not list(out.iterdir())
