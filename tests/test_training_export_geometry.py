"""export-training by run geometry (docs/TRAINING.md, step 2; docs/DESIGN_system_packs.md section 10; P5 acceptance 8).

``export_training_pairs(orientations=None)`` and the CLI default ``--orientations auto`` read each run's geometry:
a 3D run gives its xy, xz and yz slices exactly as before, and a projected run (segmentation.nuclei.project, P2a)
gives one pair of orientation ``proj``: the 2D image its segmenter saw (the recorded plane, or the max or mean
projection, read again at the run's stride, z range and timepoint) with the run's 2D labels, labels[0]. An
explicit list names the slices of 3D runs; a projected run gives its proj pair whatever the list says.

The 3D half is pinned against the base code: the manifest of the synthetic gonad (test_end_to_end's stack, run with
config/config.yaml and the classical segmenter) is the text the export wrote before this change, byte for byte
(LF line endings), and export.json has the same keys, "orientations" ["xy", "xz", "yz"] and no new key. The pinned
rows were printed by the unchanged export of commit 4059e71 (the train's base 8eb73de for segment/training.py).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile
from test_end_to_end import CHANNELS as GONAD_CHANNELS
from test_end_to_end import SPACING as GONAD_SPACING
from test_end_to_end import _synthetic_stack
from test_projected_run import IMAGE as YEAST, IMAGE_ID as YEAST_ID, Z_FOCUS, recipe, yeast_like_stack

from germquant import __version__, cli, pipeline
from germquant.config import load_config
from germquant.io.nd2_reader import Stack
from germquant.segment import training as TR

ROOT = Path(__file__).resolve().parents[1]
GONAD = "20251105_n2_nohs_HERM_001.nd2"
GONAD_ID = "20251105_n2_nohs_HERM_001"
# test_end_to_end.patched_reader's header, kept as it is so the run is the one every training test makes
GONAD_META = {"sizes": {"C": 3, "Z": 20, "Y": 64, "X": 160}, "spacing": GONAD_SPACING, "spacing_ok": True,
              "channel_names": GONAD_CHANNELS, "dtype": "float32", "is_2d": False}

# (orientation, slice_index, n_labels, z_zoom, height, width) of the base export of the synthetic gonad
BASE_DEFAULT = [("xy", 3, 28, "1.0", 80, 320), ("xy", 8, 28, "1.0", 80, 320), ("xy", 12, 28, "1.0", 80, 320),
                ("xz", 25, 28, "2.0", 40, 320), ("xz", 45, 28, "2.0", 40, 320), ("xz", 53, 28, "2.0", 40, 320),
                ("yz", 11, 1, "2.0", 40, 80), ("yz", 58, 1, "2.0", 40, 80), ("yz", 179, 1, "2.0", 40, 80)]
# the same with orientations ("yz", "xy"), slices_per_image 2 and seed 4: the given order draws first
BASE_YZ_XY = [("yz", 46, 1, "2.0", 40, 80), ("yz", 119, 1, "2.0", 40, 80), ("xy", 6, 28, "1.0", 80, 320),
              ("xy", 10, 28, "1.0", 80, 320)]
BASE_EXPORT_KEYS = ["results_dir", "orientations", "slices_per_image", "seed", "zone_column", "min_labels",
                    "max_bytes", "bytes_written", "stopped_at_cap", "runs", "germquant_version", "date"]
BASE_RUN_KEYS = ["folder", "image_id", "series_id", "source_file", "xy_stride", "z_range", "t_index",
                 "channel_role", "channel_index", "labels_path", "nuclei_path"]
HEADER = ("image_id,series_id,t_index,orientation,slice_index,n_labels,zone,image,labels,source_file,xy_stride,"
          "z_range,channel_role,channel_index,pixel_um,z_zoom,height,width")


def _base_manifest(rows) -> str:
    lines = [HEADER]
    for o, i, n, zoom, h, w in rows:
        stem = f"{GONAD_ID}_{o}_{i:04d}.tif"
        lines.append(f"{GONAD_ID},{GONAD_ID},,{o},{i},{n},,image/{stem},labels/{stem},{GONAD},1,,dna,0,0.2,{zoom},"
                     f"{h},{w}")
    return "\n".join(lines) + "\n"


def _manifest_text(pairs: Path) -> str:
    return (pairs / "manifest.csv").read_bytes().replace(b"\r\n", b"\n").decode("ascii")


def _export_json(pairs: Path) -> dict:
    return json.loads((pairs / "export.json").read_text(encoding="utf-8"))


_STACKS: dict = {}


def _sources() -> dict:
    """{file name: (header, stack)} of the two synthetic images, made once."""
    if not _STACKS:
        _STACKS[GONAD] = (GONAD_META, _synthetic_stack())
        stack, _planted = yeast_like_stack()
        c, z, y, x = stack.data.shape
        _STACKS[YEAST] = ({"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": stack.spacing, "spacing_ok": True,
                           "channel_names": list(stack.channel_names), "dtype": "float32", "is_2d": False}, stack)
    return _STACKS


def _patch_readers(mp) -> None:
    """The .nd2 reader seam, by file name, honouring xy_stride and z_range as nd2_reader.read_stack does (the
    run and the export's re-read then see the same grid)."""
    src = _sources()

    def read_meta(path, **kw):
        return src[Path(str(path)).name][0]

    def read_stack(path, *, xy_stride=1, z_range=None, **kw):
        st = src[Path(str(path)).name][1]
        data = st.data if z_range is None else st.data[:, int(z_range[0]):int(z_range[1])]
        s = int(xy_stride)
        dz, dy, dx = st.spacing
        return Stack(data=np.ascontiguousarray(data[:, :, ::s, ::s]), spacing=(dz, dy * s, dx * s),
                     channel_names=list(st.channel_names), path=Path(str(path)), spacing_ok=True)

    mp.setattr(pipeline, "read_nd2_metadata", read_meta)
    mp.setattr(pipeline, "read_stack", read_stack)


@pytest.fixture
def readers(monkeypatch):
    _patch_readers(monkeypatch)


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    """The finished runs: mixed/gonad (3D, the synthetic gonad), mixed/yeast (projected, the sharpest plane),
    max and mean (the projections), strided (the sharpest plane at xy stride 2 and z range 1 to 9)."""
    root = tmp_path_factory.mktemp("geometry")
    with pytest.MonkeyPatch.context() as mp:
        _patch_readers(mp)
        cfg = load_config(ROOT / "config" / "config.yaml")
        cfg._data["segmentation"]["nuclei"]["method"] = "classical"
        pipeline.process_image(GONAD, cfg, root / "mixed" / "gonad", prov=None)
        pipeline.process_image(YEAST, recipe(**{"spots.enabled": False}), root / "mixed" / "yeast", prov=None)
        for mode in ("max", "mean"):
            pipeline.process_image(YEAST, recipe(**{"spots.enabled": False, "segmentation.nuclei.project": mode}),
                                   root / mode, prov=None)
        pipeline.process_image(YEAST, recipe(**{"spots.enabled": False}), root / "strided", prov=None,
                               xy_stride=2, z_range=(1, 9))
    return root


def _labels(folder: Path, iid: str) -> np.ndarray:
    return tifffile.imread(folder / f"{iid}__nuclei_labels.tif")


def _u16(a) -> np.ndarray:
    return np.clip(np.rint(np.asarray(a, dtype=np.float64)), 0, 65535).astype(np.uint16)


# ---- 3D runs: unchanged -----------------------------------------------------------------------------------

def test_3d_export_manifest_and_export_json_unchanged(tmp_path, results, readers):
    gonad = results / "mixed" / "gonad"
    assert TR.ORIENTATIONS == ("xy", "xz", "yz")              # the plugin's default export list imports it
    auto = TR.export_training_pairs(gonad, tmp_path / "auto")
    assert _manifest_text(tmp_path / "auto") == _base_manifest(BASE_DEFAULT)
    assert list(auto.columns) == TR.MANIFEST_COLS and len(auto) == 9
    exp = _export_json(tmp_path / "auto")
    assert list(exp) == BASE_EXPORT_KEYS                       # no new key for a 3D-only tree
    assert exp["orientations"] == ["xy", "xz", "yz"]
    assert (exp["slices_per_image"], exp["seed"], exp["zone_column"], exp["min_labels"]) == (3, 0, "zone", 3)
    assert exp["max_bytes"] == 2e9 and exp["stopped_at_cap"] is False and exp["germquant_version"] == __version__
    assert isinstance(exp["bytes_written"], int) and exp["bytes_written"] > 0
    assert exp["results_dir"] == str(gonad) and len(exp["runs"]) == 1
    run = exp["runs"][0]
    assert list(run) == BASE_RUN_KEYS
    assert run == {"folder": str(gonad), "image_id": GONAD_ID, "series_id": GONAD_ID, "source_file": GONAD,
                   "xy_stride": 1, "z_range": None, "t_index": None, "channel_role": "dna", "channel_index": 0,
                   "labels_path": str(gonad / f"{GONAD_ID}__nuclei_labels.tif"),
                   "nuclei_path": str(gonad / f"{GONAD_ID}__nuclei.csv")}
    # the old default passed explicitly, "auto" and the CLI default write the same pairs
    TR.export_training_pairs(gonad, tmp_path / "explicit", orientations=("xy", "xz", "yz"))
    TR.export_training_pairs(gonad, tmp_path / "auto_text", orientations="auto")
    assert cli.main(["export-training", str(gonad), "--out", str(tmp_path / "cli")]) == 0
    assert cli.main(["export-training", str(gonad), "--out", str(tmp_path / "cli_list"),
                     "--orientations", "xy,xz,yz"]) == 0
    for name in ("explicit", "auto_text", "cli", "cli_list"):
        assert _manifest_text(tmp_path / name) == _base_manifest(BASE_DEFAULT), name
        e = _export_json(tmp_path / name)
        assert list(e) == BASE_EXPORT_KEYS and e["orientations"] == ["xy", "xz", "yz"], name
    for r in auto.to_dict("records"):
        for name in ("explicit", "cli"):
            for kind in ("image", "labels"):
                assert np.array_equal(tifffile.imread(tmp_path / "auto" / r[kind]),
                                      tifffile.imread(tmp_path / name / r[kind]))
    # an explicit list keeps its order (it decides the seeded draw) in the slices and in export.json
    TR.export_training_pairs(gonad, tmp_path / "yz_xy", orientations=("yz", "xy"), slices_per_image=2, seed=4)
    assert _manifest_text(tmp_path / "yz_xy") == _base_manifest(BASE_YZ_XY)
    assert _export_json(tmp_path / "yz_xy")["orientations"] == ["yz", "xy"]
    # the pixels: the raw dna plane as uint16 and the run's own labels
    stack, labels = _synthetic_stack(), _labels(gonad, GONAD_ID)
    for r in auto[auto["orientation"] == "xy"].to_dict("records"):
        k = int(r["slice_index"])
        assert np.array_equal(tifffile.imread(tmp_path / "auto" / r["image"]), _u16(stack.data[0][k]))
        assert np.array_equal(tifffile.imread(tmp_path / "auto" / r["labels"]), labels[k])
    for bad in (("zz",), ("proj",), ("xy", "proj"), (), "xy,bogus"):
        with pytest.raises(ValueError, match="orientations"):
            TR.export_training_pairs(gonad, tmp_path / "bad", orientations=bad)


# ---- projected runs: the 2D image the segmenter saw ----------------------------------------------------------

def test_projected_run_exports_proj_pairs(tmp_path, results, readers):
    yeast = results / "mixed" / "yeast"
    stack = _sources()[YEAST][1]
    labels = _labels(yeast, YEAST_ID)
    record = {"mode": "plane", "plane": "sharpest", "z": Z_FOCUS}
    assert json.loads((yeast / f"{YEAST_ID}__stages.json").read_text(encoding="utf-8"))["segment_projection"] == record
    man = TR.export_training_pairs(yeast, tmp_path / "plane")
    assert list(man.columns) == TR.MANIFEST_COLS and len(man) == 1
    r = man.iloc[0].to_dict()
    stem = f"{YEAST_ID}_proj_{Z_FOCUS:04d}.tif"
    assert (r["orientation"], int(r["slice_index"]), r["image"], r["labels"]) == ("proj", Z_FOCUS, f"image/{stem}",
                                                                                f"labels/{stem}")
    assert (r["channel_role"], int(r["channel_index"]), float(r["pixel_um"]), float(r["z_zoom"])) == \
        ("transmitted", 0, 0.1, 1.0)
    assert (int(r["height"]), int(r["width"]), int(r["n_labels"])) == (110, 200, 8)
    im = tifffile.imread(tmp_path / "plane" / r["image"])
    lab = tifffile.imread(tmp_path / "plane" / r["labels"])
    assert im.dtype == np.uint16 and lab.dtype == np.int32
    assert np.array_equal(im, _u16(stack.data[0][Z_FOCUS]))           # the plane the run segmented
    assert np.array_equal(lab, labels[0])                               # the run's 2D labels
    exp = _export_json(tmp_path / "plane")
    assert exp["orientations"] == ["proj"]
    assert exp["run_orientations"] == [{"image_id": YEAST_ID, "geometry": "projected", "orientations": ["proj"],
                                        "segment_projection": record}]
    assert list(exp) == BASE_EXPORT_KEYS[:2] + ["run_orientations"] + BASE_EXPORT_KEYS[2:]
    # an explicit list names the slices of 3D runs: a projected run still gives its proj pair
    again = TR.export_training_pairs(yeast, tmp_path / "listed", orientations=("xy", "xz", "yz"))
    pd.testing.assert_frame_equal(again, man)
    assert _manifest_text(tmp_path / "listed") == _manifest_text(tmp_path / "plane")
    assert _export_json(tmp_path / "listed")["orientations"] == ["proj"]
    # the projections: one image, index 0
    for mode in ("max", "mean"):
        m = TR.export_training_pairs(results / mode, tmp_path / mode)
        assert len(m) == 1 and m["orientation"].iloc[0] == "proj" and int(m["slice_index"].iloc[0]) == 0
        img = stack.data[0]
        want = img.max(axis=0) if mode == "max" else img.mean(axis=0, dtype=np.float64).astype(np.float32)
        assert np.array_equal(tifffile.imread(tmp_path / mode / m["image"].iloc[0]), _u16(want))
        assert np.array_equal(tifffile.imread(tmp_path / mode / m["labels"].iloc[0]),
                              _labels(results / mode, YEAST_ID)[0])
        rec = _export_json(tmp_path / mode)["run_orientations"][0]["segment_projection"]
        assert rec == {"mode": mode, "plane": None, "z": None}
    # read again at the run's geometry: xy stride 2 and z range 1 to 9, the plane counted inside that range
    s = TR.export_training_pairs(results / "strided", tmp_path / "strided")
    rec = _export_json(tmp_path / "strided")["run_orientations"][0]["segment_projection"]
    assert rec["z"] == Z_FOCUS - 1 and int(s["slice_index"].iloc[0]) == Z_FOCUS - 1
    assert (int(s["height"].iloc[0]), int(s["width"].iloc[0]), float(s["pixel_um"].iloc[0])) == (55, 100, 0.2)
    exp_s = _export_json(tmp_path / "strided")["runs"][0]
    assert exp_s["xy_stride"] == 2 and exp_s["z_range"] == [1, 9]
    assert np.array_equal(tifffile.imread(tmp_path / "strided" / s["image"].iloc[0]),
                          _u16(stack.data[0][1:9, ::2, ::2][Z_FOCUS - 1]))
    # the pairs load back with their orientation, with or without the manifest
    pairs = TR.load_pairs(tmp_path / "plane")
    assert [(p["image_id"], p["orientation"], p["slice_index"]) for p in pairs] == [(YEAST_ID, "proj", Z_FOCUS)]
    (tmp_path / "plane" / "manifest.csv").unlink()
    by_name = TR.load_pairs(tmp_path / "plane")
    assert [(p["image_id"], p["orientation"], p["slice_index"]) for p in by_name] == [(YEAST_ID, "proj", Z_FOCUS)]


def test_a_mixed_tree_exports_each_run_by_its_geometry(tmp_path, results, readers):
    mixed = results / "mixed"
    man = TR.export_training_pairs(mixed, tmp_path / "pairs")
    assert sorted(set(man["orientation"])) == ["proj", "xy", "xz", "yz"]
    exp = _export_json(tmp_path / "pairs")
    assert exp["orientations"] == ["xy", "xz", "yz", "proj"]          # the union, slices first
    assert [r["image_id"] for r in exp["runs"]] == [GONAD_ID, YEAST_ID]
    assert exp["run_orientations"] == [
        {"image_id": GONAD_ID, "geometry": "3d", "orientations": ["xy", "xz", "yz"]},
        {"image_id": YEAST_ID, "geometry": "projected", "orientations": ["proj"],
         "segment_projection": {"mode": "plane", "plane": "sharpest", "z": Z_FOCUS}}]
    # the 3D run's slices do not depend on the other runs in the tree (the draw is seeded per image)
    TR.export_training_pairs(mixed / "gonad", tmp_path / "alone")
    lines = _manifest_text(tmp_path / "pairs").splitlines()
    assert [ln for ln in lines if ln.startswith(GONAD_ID + ",")] == _manifest_text(tmp_path / "alone").splitlines()[1:]
    assert len(man) == 10 and len([ln for ln in lines if ln.startswith(YEAST_ID + ",")]) == 1
    # an explicit 3D subset: the slices of the 3D run, the proj pair of the projected one
    sub = TR.export_training_pairs(mixed, tmp_path / "xy_only", orientations=("xy",))
    assert list(sub["orientation"]) == ["xy", "xy", "xy", "proj"]
    assert _export_json(tmp_path / "xy_only")["orientations"] == ["xy", "proj"]
    # the command line: auto by default, a list on request
    assert cli.main(["export-training", str(mixed), "--out", str(tmp_path / "cli")]) == 0
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / "cli" / "manifest.csv"),
                                  pd.read_csv(tmp_path / "pairs" / "manifest.csv"))
    assert cli.main(["export-training", str(mixed), "--out", str(tmp_path / "cli_yz"), "--orientations", "yz"]) == 0
    assert list(pd.read_csv(tmp_path / "cli_yz" / "manifest.csv")["orientation"]) == ["yz", "yz", "yz", "proj"]
