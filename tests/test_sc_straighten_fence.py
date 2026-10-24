"""The fence around sc-straighten: its modules import no napari, Qt, torch or CuPy; `germquant collect` on a tree
without sc_straighten records writes exactly the old file set (and stacks the new tables when records exist); the
napari layer helpers return valid LayerData tuples without napari."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_sc_straighten_cli import _copy, _run, finished, finished_run  # noqa: F401  (fixture)

from germquant import batch
from germquant.sc.straighten import layers, outputs
from germquant.sc.straighten import schema as S

ROOT = Path(__file__).resolve().parents[1]


def test_imports_stay_light():
    code = ("import sys; import germquant.sc.straighten.runner, germquant.sc.straighten.layers, "
            "germquant.sc.straighten.pipeline, germquant.cli_sc_straighten, germquant.benchmark.sc_ribbons; "
            "bad = sorted(m for m in sys.modules if m.split('.')[0] in "
            "('napari', 'PyQt5', 'PyQt6', 'PySide2', 'PySide6', 'qtpy', 'torch', 'cupy', 'cellpose', 'spotmax')); "
            "print(','.join(bad))")
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, check=True)
    assert out.stdout.strip() == "", out.stdout


@pytest.fixture(scope="module")
def finished_plain(tmp_path_factory):
    """The same gonad under an ordinary image id (no trailing underscore), which `germquant collect` stacks."""
    return finished(tmp_path_factory, "20251105_N2_nohs_HERM_001.nd2")


def test_collect_writes_the_old_file_set_without_records(finished_run, finished_plain, tmp_path):  # noqa: F811
    for fr, sub in ((finished_run, "trailing"), (finished_plain, "plain")):
        root = _copy(fr, tmp_path / sub)
        before = {p.name for p in root.iterdir()}
        counts = batch.collect(root)
        assert not any(k.startswith("sc_") for k in counts if k not in batch.schema.TABLES)
        assert not [n for n in S.BATCH_NAMES.values() if (root / n).exists()]
        new = {p.name for p in root.iterdir()} - before
        assert new and all(n.startswith("batch_") for n in new)
    # the ordinary id is stacked into the main tables (an id ending in an underscore is not: a known issue of
    # batch.collect outside this stage)
    assert counts["nuclei"] > 0 and len(pd.read_csv(root / "batch_summary.csv")) == 1


def test_collect_skips_stale_records_and_follows_image_ids(finished_plain, tmp_path):
    _r, iid, stack = finished_plain
    root = _copy(finished_plain, tmp_path)
    assert _run(root, stack)[0].status == "done"
    assert outputs.collect_tables(root, image_ids={"some other image"}) is None
    c = outputs.collect_tables(root, image_ids={iid})
    assert c["sc_strands"] > 0
    nuc = root / "img" / f"{iid}__nuclei.csv"
    nuc.write_bytes(nuc.read_bytes() + b"\n")              # the germquant run changed after sc-straighten
    c = outputs.collect_tables(root)
    assert c["sc_strands"] == 0
    summ = pd.read_csv(root / S.BATCH_NAMES["summary"])
    assert list(summ["status"]) == ["stale"]


def test_collect_stacks_the_new_tables_with_records(finished_run, tmp_path):  # noqa: F811
    _r, iid, stack = finished_run
    root = _copy(finished_run, tmp_path)
    _run(root, stack)
    counts = batch.collect(root)
    for name in S.BATCH_NAMES.values():
        assert (root / name).is_file(), name
    assert counts["sc_strands"] > 0 and counts["sc_straighten_summary"] == 1


def test_layers_without_napari(finished_run, tmp_path):  # noqa: F811
    _r, iid, stack = finished_run
    root = _copy(finished_run, tmp_path)
    _run(root, stack)
    folder = root / "img"
    lay = layers.strand_layers(folder, iid)
    assert lay and all(len(t) == 3 and t[2] in ("shapes", "points", "image") for t in lay)
    paths = lay[0][0]
    assert paths and all(p.ndim == 2 and p.shape[1] == 3 for p in paths)
    import pandas as pd

    nid = int(pd.read_csv(folder / f"{iid}{S.SUFFIXES['strands']}")["nucleus_id"].iloc[0])
    gal, strands, desc = layers.strip_gallery(folder, iid, nid)
    assert gal.ndim == 3 and gal.shape[0] == len(strands) and desc["format"] == "germquant.sc_straighten.strips"
    sl = layers.strip_layers(folder, iid, nid)
    assert sl[0][2] == "image" and np.isfinite(sl[0][0]).any()
    il = layers.image_layers(np.zeros((3, 4, 5)), (1, 2, 3), (0.2, 0.1, 0.1))
    assert np.allclose(il[0][1]["translate"], (0.2, 0.2, 0.3)) and il[0][2] == "image"
