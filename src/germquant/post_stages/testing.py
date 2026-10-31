"""The fixture post-stage that the engine and plugin tests share (plugin need I14; docs/POST_STAGES.md).

``fixture_nuclei`` measures, for each nucleus of the run's label image, its voxel count and volume, a statistic of
the foci channel inside it (mean or max) times a factor, and its spot count from the spots table. It reads every
kind of input (a channel role, two tables, the label image, the image) and declares one parameter of each type
the parameter form shows (float, int, choice, bool); ``fail: true`` makes it raise, for the failure tests.

Both forms come from the same source text below, so nothing is copied and nothing can drift:
* ``register_fixture()`` registers it in this process (the engine's in-process tests);
* ``write_fixture_distribution(folder)`` writes a dist-info with the entry point plus the package into ``folder``,
  which any interpreter that has ``folder`` on sys.path (PYTHONPATH) then discovers from metadata: the case of a
  subprocess engine, which cannot see an in-process registration.

``write_distribution`` writes any such distribution, for tests of broken declarations. Nothing here imports more
than the standard library and yaml until the fixture runs.
"""
from __future__ import annotations

import re
import types
from pathlib import Path

FIXTURE_NAME = "fixture_nuclei"
FIXTURE_DIST = "germquant-fixture-post-stage"
FIXTURE_PACKAGE = "germquant_fixture_post_stage"
FIXTURE_VERSION = "1.0.0"
FIXTURE_TABLE = "per_nucleus"
FIXTURE_COLUMNS = ("nucleus_id", "n_voxels", "volume_um3", "foci_value", "foci_value_scaled", "n_spots")

_DECLARATION = """format: 1
post_stages:
  {name}:
    version: "1.0.0"
    title: "Fixture: per-nucleus voxels, volume, foci statistic and spot count"
    reads: {{roles: [foci], tables: [nuclei, spots], labels: true, image: true}}
    writes:
      per_nucleus: [nucleus_id, n_voxels, volume_um3, foci_value, foci_value_scaled, n_spots]
    params:
      - {{key: factor, type: float, min: 0.1, max: 10.0, default: 1.0, label: "Scale factor of foci_value"}}
      - {{key: min_voxels, type: int, min: 0, max: 1000000000, default: 0, label: "Smallest nucleus (voxels)"}}
      - {{key: statistic, type: choice, choices: [mean, max], default: mean, label: "Foci statistic"}}
      - {{key: fail, type: bool, default: false, label: "Fail on purpose (tests)"}}
"""

_SOURCE = '''"""germquant's fixture post-stage (germquant.post_stages.testing, docs/POST_STAGES.md).

Its declaration is germquant_post_stage.yaml next to this file. numpy and pandas are imported when it runs,
never when the package is imported.
"""

COLUMNS = ["nucleus_id", "n_voxels", "volume_um3", "foci_value", "foci_value_scaled", "n_spots"]


def measure(inputs, params):
    """Per nucleus: voxels, volume, the foci mean or max (times params["factor"]) and the spot count."""
    import numpy as np
    import pandas as pd

    if params["fail"]:
        raise RuntimeError("the fixture post-stage was asked to fail (fail: true)")
    labels = np.asarray(inputs.labels(), dtype=np.int64)
    foci = np.asarray(inputs.image("foci"), dtype=np.float64)
    if foci.shape != labels.shape:
        raise ValueError(f"the foci channel {foci.shape} and the labels {labels.shape} differ in shape")
    nuclei = inputs.tables("nuclei")
    spots = inputs.tables("spots")
    ids = sorted(int(v) for v in nuclei["nucleus_id"])
    flat = labels.ravel()
    values = foci.ravel()
    size = max([0, *ids, int(flat.max()) if flat.size else 0]) + 1
    counts = np.bincount(flat, minlength=size)
    sums = np.bincount(flat, weights=values, minlength=size)
    maxes = np.full(size, -np.inf)
    np.maximum.at(maxes, flat, values)
    n_spots = spots["nucleus_id"].astype(int).value_counts().to_dict() if len(spots) else {}
    spacing = inputs.spacing
    voxel = float(spacing[0] * spacing[1] * spacing[2]) if spacing else float("nan")
    rows = []
    for i in ids:
        n = int(counts[i])
        if n < params["min_voxels"]:
            continue
        if n == 0:
            value = float("nan")
        elif params["statistic"] == "max":
            value = float(maxes[i])
        else:
            value = float(sums[i] / n)
        rows.append({"nucleus_id": i, "n_voxels": n, "volume_um3": n * voxel, "foci_value": value,
                     "foci_value_scaled": value * params["factor"], "n_spots": int(n_spots.get(i, 0))})
    return {"per_nucleus": pd.DataFrame(rows, columns=COLUMNS)}
'''


def declaration_text(name: str = FIXTURE_NAME) -> str:
    """The fixture's germquant_post_stage.yaml, declaring it under ``name`` (other names serve tests)."""
    return _DECLARATION.format(name=name)


def fixture_declaration() -> dict:
    """The fixture's declaration (the mapping under post_stages.fixture_nuclei), a fresh copy each call."""
    import yaml

    return yaml.safe_load(declaration_text())["post_stages"][FIXTURE_NAME]


def fixture_function():
    """The fixture's callable, made from the same source text the on-disk package holds (never put in
    sys.modules)."""
    module = types.ModuleType(f"{FIXTURE_PACKAGE}_in_process")
    exec(compile(_SOURCE, f"<{FIXTURE_PACKAGE}>", "exec"), module.__dict__)  # noqa: S102 - the engine's own text
    return module.measure


def register_fixture(*, replace: bool = False) -> str:
    """Register the fixture in this process (germquant.post_stages.register_post_stage); returns its name."""
    from .discovery import register_post_stage

    register_post_stage(FIXTURE_NAME, fixture_function(), fixture_declaration(), replace=replace)
    return FIXTURE_NAME


def write_distribution(folder, *, name: str, version: str, entry_points: dict[str, str],
                       files: dict[str, str | bytes]) -> Path:
    """Write an installed distribution into ``folder`` (a sys.path entry): <name>-<version>.dist-info with
    METADATA, top_level.txt and the ``entry_points`` {post-stage name: "module:attr"} in the germquant.post_stages
    group, plus ``files`` {relative path: text or bytes}. Returns ``folder``."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    info = folder / f"{re.sub(r'[-_.]+', '_', name)}-{version}.dist-info"
    info.mkdir(exist_ok=True)
    (info / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n"
                                   f"Summary: a germquant post-stage written by germquant.post_stages.testing\n",
                                   encoding="utf-8", newline="\n")
    lines = "".join(f"{k} = {v}\n" for k, v in entry_points.items())
    (info / "entry_points.txt").write_text(f"[germquant.post_stages]\n{lines}", encoding="utf-8", newline="\n")
    tops = sorted({rel.split("/", 1)[0].removesuffix(".py") for rel in files})
    (info / "top_level.txt").write_text("".join(f"{t}\n" for t in tops), encoding="utf-8", newline="\n")
    (info / "RECORD").write_text("", encoding="utf-8")
    for rel, content in files.items():
        path = folder / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(content, bytes):
            path.write_bytes(content)
        else:
            path.write_text(content, encoding="utf-8", newline="\n")
    return folder


def write_fixture_distribution(folder) -> Path:
    """The fixture as an on-disk distribution in ``folder``: germquant_fixture_post_stage-1.0.0.dist-info (entry
    point fixture_nuclei = germquant_fixture_post_stage:measure) and the package germquant_fixture_post_stage/
    with its declaration file. Put ``folder`` on PYTHONPATH (or sys.path) to install it; returns ``folder``."""
    return write_distribution(folder, name=FIXTURE_DIST, version=FIXTURE_VERSION,
                              entry_points={FIXTURE_NAME: f"{FIXTURE_PACKAGE}:measure"},
                              files={f"{FIXTURE_PACKAGE}/__init__.py": _SOURCE,
                                     f"{FIXTURE_PACKAGE}/germquant_post_stage.yaml": declaration_text()})
