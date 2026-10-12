"""Pack goldens, recipe tier (docs/DESIGN_system_packs.md, section 6.9; gate G8).

Each recipe of a pack runs on the pack's synthetic example with segmentation.nuclei.method forced to classical
by mutating ``cfg._data`` (as tests/test_golden.py does), so the config_hash sentinel is the recipe's true
hash, and its tables are compared cell for cell with ``example/expected/recipe/<recipe>/`` (nuclei.csv,
spots.csv, image_summary.csv, config_hash.txt). Classical runs never resolve a pack model, so the tier runs
with an empty model store (``$GERMQUANT_MODELS`` pointed at an empty folder for the run). CI has no bioio, so
``read_example`` reads the OME-TIFF with tifffile and process_image takes it through its ``read_fn`` seam;
where bioio exists, a test shows both readers give the same Stack.

The model tier (a fetched public crop on the RTX 5090, compared exactly, twice) arrives with the 1.0.0 packs.
This module imports the pipeline, and so torch; germquant.packs does not import it.
"""
from __future__ import annotations

import contextlib
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .io.ome_header import ome_fields as _ome_fields

TABLES = ("nuclei", "spots", "image_summary")
VOLATILE = ["git_sha", "run_timestamp", "pipeline_version", "file_path", "config_hash"]


def _to_czyx(data: np.ndarray, axes: str, t: int | None) -> np.ndarray:
    axes = axes.upper()
    for ax in list(axes):
        if ax not in "TCZYX":
            if data.shape[axes.index(ax)] != 1:
                raise ValueError(f"example axis {ax} of size {data.shape[axes.index(ax)]} is not supported")
            data = np.take(data, 0, axis=axes.index(ax))
            axes = axes.replace(ax, "")
    if "T" in axes:
        n_t = data.shape[axes.index("T")]
        if t is None and n_t > 1:
            raise ValueError(f"the example has {n_t} time points; pass t")
        data = np.take(data, int(t or 0), axis=axes.index("T"))
        axes = axes.replace("T", "")
    for ax in "CZ":
        if ax not in axes:
            data = data[np.newaxis]
            axes = ax + axes
    order = [axes.index(a) for a in "CZYX"]
    return np.ascontiguousarray(np.transpose(data, order))


def read_example(path, *, xy_stride: int = 1, z_range=None, t: int | None = None):
    """``(meta, Stack)`` of an OME-TIFF example read with tifffile, the contract of the bioio branch of
    pipeline._read_input (the same spacing rule: a missing or non-positive size gives (1, 1, 1) with
    spacing_ok False)."""
    import tifffile

    from .io.nd2_reader import Stack

    path = Path(path)
    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        axes = series.axes
        raw = series.asarray()
        names, (pz, py, px) = _ome_fields(tf.ome_metadata)
    data = _to_czyx(raw, axes, t)
    c, z, y, x = data.shape
    spacing_ok = all(v is not None for v in (pz, py, px))
    file_spacing = (float(pz), float(py), float(px)) if spacing_ok else (1.0, 1.0, 1.0)
    dz, dy, dx = file_spacing
    if z_range is not None:
        data = data[:, z_range[0]:z_range[1]]
    if xy_stride > 1:
        data = data[:, :, ::xy_stride, ::xy_stride]
        dy, dx = dy * xy_stride, dx * xy_stride
    data = np.ascontiguousarray(data)
    if len(names) != c:
        names = (names + [str(i) for i in range(len(names), c)])[:c]
    meta = {"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": file_spacing, "spacing_ok": spacing_ok,
            "channel_names": list(names), "dtype": str(raw.dtype), "is_2d": z == 1, "n_timepoints": 1,
            "period_s": None, "t_seconds": None, "n_scenes": 1, "reader": "tifffile"}
    stack = Stack(data=data, spacing=(dz, dy, dx), channel_names=list(names), path=path,
                  downsample_xy=xy_stride, spacing_ok=spacing_ok, t_index=t)
    return meta, stack


@contextlib.contextmanager
def _empty_models_dir():
    """``$GERMQUANT_MODELS`` pointed at a new empty folder for the duration (restored afterwards)."""
    old = os.environ.get("GERMQUANT_MODELS")
    with tempfile.TemporaryDirectory(prefix="gq-empty-models-") as tmp:
        os.environ["GERMQUANT_MODELS"] = tmp
        try:
            yield Path(tmp)
        finally:
            if old is None:
                os.environ.pop("GERMQUANT_MODELS", None)
            else:
                os.environ["GERMQUANT_MODELS"] = old


def example_image(pack, recipe) -> Path:
    """The example a recipe golden runs on: the recipe's ``example`` id, else the pack's first local example."""
    want = (recipe.raw or {}).get("example")
    for ex in pack.data.get("examples") or []:
        if not isinstance(ex, dict) or not ex.get("image"):
            continue
        if want is None or ex.get("id") == want:
            return pack.folder / str(ex["image"])
    raise FileNotFoundError(f"pack {pack.ref}: no local example image for recipe {recipe.id}")


def expected_dir(pack, recipe_id: str) -> Path:
    return pack.folder / "example" / "expected" / "recipe" / recipe_id


def _read(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, low_memory=False, float_precision="round_trip")
    return df.drop(columns=[c for c in VOLATILE if c in df.columns])


def _tolerance() -> dict:
    rtol = os.environ.get("GERMQUANT_GOLDEN_RTOL")
    return {"check_exact": False, "rtol": float(rtol), "atol": 0.0} if rtol else {"check_exact": True}


def compare_table(got: pd.DataFrame, ref: pd.DataFrame, name: str) -> list[str]:
    """The comparison of tests/test_golden.py: the reference columns are all present in the same relative
    order (added columns are allowed), the row count matches, numbers are equal (GERMQUANT_GOLDEN_RTOL in CI)
    and strings are identical."""
    problems = []
    common = [c for c in got.columns if c in set(ref.columns)]
    if common != list(ref.columns):
        return [f"{name}: columns renamed, removed or reordered"]
    if len(got) != len(ref):
        return [f"{name}: {len(got)} rows, expected {len(ref)}"]
    tol = _tolerance()
    for col in ref.columns:
        a, b = got[col], ref[col]
        if pd.api.types.is_numeric_dtype(b):
            try:
                pd.testing.assert_series_equal(a.astype(float), b.astype(float), check_names=False, **tol)
            except AssertionError:
                problems.append(f"{name}.{col} differs")
        elif not (a.astype("string").fillna("<NA>") == b.astype("string").fillna("<NA>")).all():
            problems.append(f"{name}.{col} differs")
    return problems


def run_recipe_golden(pack, recipe: str | None = None, *, update: bool = False, work_dir=None) -> dict:
    """Run one recipe golden (module docstring). ``update`` writes the expected tables and the sentinel
    instead of comparing. Returns ``{"recipe", "ok", "problems", "config_hash", "out"}``."""
    from .config import load_config
    from .packs import recipe_profile_path
    from .pipeline import process_image

    r = pack.recipe(recipe)
    problems: list[str] = []
    cfg = load_config(recipe_profile_path(pack, r))
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"     # not an override: the hash stays the recipe's
    if cfg.hash != r.expected_config_hash:
        problems.append(f"config_hash {cfg.hash} != the pack's expected_config_hash {r.expected_config_hash}")
    image = example_image(pack, r)
    gold = expected_dir(pack, r.id)
    with contextlib.ExitStack() as stack:
        work = Path(work_dir) if work_dir is not None else Path(stack.enter_context(tempfile.TemporaryDirectory()))
        stack.enter_context(_empty_models_dir())
        out = work / f"golden_{pack.id}_{r.id}"
        res = process_image(image, cfg, out, prov=None, read_fn=read_example)
        files = {t: out / f"{res['image_id']}__{t}.csv" for t in TABLES}
        if update:
            gold.mkdir(parents=True, exist_ok=True)
            for t, f in files.items():
                _read(f).to_csv(gold / f"{t}.csv", index=False, lineterminator="\n")
            (gold / "config_hash.txt").write_bytes((cfg.hash + "\n").encode("ascii"))
            return {"recipe": r.id, "ok": not problems, "problems": problems, "config_hash": cfg.hash,
                    "out": str(out)}
        if not gold.is_dir():
            return {"recipe": r.id, "ok": False, "problems": [f"no expected folder {gold}"], "config_hash": cfg.hash,
                    "out": str(out)}
        sentinel = (gold / "config_hash.txt").read_text(encoding="ascii").strip()
        if sentinel != cfg.hash:
            problems.append(f"sentinel config_hash.txt {sentinel} != the recipe's hash {cfg.hash}")
        for t, f in files.items():
            ref = pd.read_csv(gold / f"{t}.csv", low_memory=False, float_precision="round_trip")
            problems += compare_table(_read(f), ref, t)
    return {"recipe": r.id, "ok": not problems, "problems": problems, "config_hash": cfg.hash, "out": str(out)}


def run_pack_goldens(pack, recipe: str | None = None, *, update: bool = False) -> list[dict]:
    """Every recipe golden of a pack, or the one ``recipe`` names. With ``update``, every recipe with an example
    is written. Without it, every recipe of an outside pack (root: pack) is compared, and one without an
    expected folder is a failed result (rule R7: every outside recipe ships a recipe golden), never skipped; a
    root: checkout pack (the C. elegans pack, whose recipes are the frozen engine files and whose goldens are
    tests/golden) runs the recipes that have an expected folder."""
    if recipe is not None:
        return [run_recipe_golden(pack, recipe, update=update)]
    out = []
    for r in pack.recipes:
        if update:
            try:
                example_image(pack, r)
            except FileNotFoundError:
                continue
        elif not expected_dir(pack, r.id).is_dir():
            if pack.root == "pack":
                gold = expected_dir(pack, r.id)
                out.append({"recipe": r.id, "ok": False, "config_hash": None, "out": None,
                            "problems": [f"no recipe golden: {gold} is missing (every recipe of an outside pack "
                                         f"ships one, rule R7; write it with `germquant pack golden {pack.id} "
                                         f"--recipe {r.id} --update`)"]})
            continue
        out.append(run_recipe_golden(pack, r.id, update=update))
    return out
