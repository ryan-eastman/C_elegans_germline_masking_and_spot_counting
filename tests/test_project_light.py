"""The project API's frozen signatures and its light entry (PLAN.md section 9.1, decision D28; gate G6, rule R5).

* germquant.project exports every name of section 9.1 with exactly the frozen parameters, so the units and
  the plugin that code against it (P5-RUN, P5-TUNE, P6, the P7 draft) can rely on them.
* In a fresh interpreter, ``import germquant.project`` loads none of torch, cupy, spotmax, cellpose, bioio
  and dask.array.
"""
from __future__ import annotations

import inspect
from pathlib import Path

from project_helpers import heavy_after

import germquant.project as P
from germquant.project import headers

_EMPTY = inspect.Parameter.empty
POS, KW = inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY

# name: [(parameter, kind, default)], from PLAN.md section 9.1 (frozen on day 1)
FROZEN = {
    "new_project": [("root", POS, _EMPTY), ("pack", POS, _EMPTY), ("recipe", POS, None), ("images", POS, ()),
                    ("name", KW, None), ("plugin_version", KW, None)],
    "open_project": [("root", POS, _EMPTY)],
    "save_project": [("project", POS, _EMPTY), ("plugin_version", KW, None), ("device", KW, None),
                     ("interpreters", KW, None), ("runner", KW, None), ("xy_stride", KW, None), ("z_range", KW, None),
                     ("overrides", KW, None), ("samples", KW, None)],
    "add_images": [("project", POS, _EMPTY), ("paths", POS, _EMPTY), ("plugin_version", KW, None)],
    "project_config": [("project", POS, _EMPTY)],
    "preview_roles": [("project", POS, _EMPTY), ("roles", POS, _EMPTY), ("image_id", KW, None)],
    "write_recipe": [("project", POS, _EMPTY), ("edits", KW, None), ("unset", KW, ()), ("roles", KW, None),
                     ("origin", KW, "form"), ("plugin_version", KW, None)],
    "revert_recipe": [("project", POS, _EMPTY), ("plugin_version", KW, None)],
    "use_model": [("project", POS, _EMPTY), ("model", POS, _EMPTY), ("card", KW, None), ("plugin_version", KW, None)],
    "run_project": [("project", POS, _EMPTY), ("images", POS, None), ("resume", POS, True), ("device", POS, None),
                    ("shard", POS, None), ("progress", POS, None), ("via", KW, "cli"), ("run_id", KW, None),
                    ("should_stop", KW, None)],
    "project_status": [("project", POS, _EMPTY)],
    "collect_project": [("project", POS, _EMPTY)],
    "export_commands": [("project", POS, _EMPTY), ("shell", KW, "sh"), ("shards", KW, None)],
    "break_lock": [("project", POS, _EMPTY), ("host", KW, _EMPTY), ("pid", KW, _EMPTY)],
    "list_rounds": [("project", POS, _EMPTY)],
    "new_round": [("project", POS, _EMPTY), ("name", POS, None), ("plugin_version", KW, None)],
}


def _params(fn) -> list[tuple]:
    return [(p.name, p.kind, p.default) for p in inspect.signature(fn).parameters.values()]


def test_the_api_has_the_frozen_signatures():
    assert set(FROZEN) | {"Project"} <= set(P.__all__)
    for name, want in FROZEN.items():
        assert _params(getattr(P, name)) == want, name
    assert _params(headers.read_header) == [("path", POS, _EMPTY)]
    fields = P.Project.__dataclass_fields__
    assert {"root", "manifest", "warnings"} <= set(fields)
    proj = P.Project(root=Path("x"), manifest={"images": [{"id": "a"}], "pack": {"id": "p"}})
    assert proj.images == [{"id": "a"}] and proj.pack == {"id": "p"} and proj.warnings == []


def test_import_germquant_project_is_light():
    heavy, r = heavy_after("import germquant.project\nimport germquant.project.headers\n")
    assert r.returncode == 0, r.stderr[-2000:]
    assert heavy == "", f"import germquant.project loaded {heavy}"
