"""export-training by run geometry (docs/TRAINING.md, step 2; docs/DESIGN_system_packs.md section 10; P5 acceptance 8).

``export_training_pairs(orientations=None)`` and the CLI default ``--orientations auto`` read each run's geometry:
a 3D run gives its xy, xz and yz slices exactly as before, and a projected run (segmentation.nuclei.project, P2a)
gives one pair of orientation ``proj``: the 2D image its segmenter saw (the recorded plane, or the max or mean
projection, read again at the run's stride, z range and timepoint) with the run's 2D labels, labels[0]. A
single-plane run (io.planar, P2b) gives its plane as xy pairs: the whole plane, or above 1024 px the fewest tiles
of at most 1024 px a side, with each tile's offset in the manifest (tile_y0, tile_x0). An explicit list names the
slices of 3D runs; a projected or single-plane run gives its 2D pairs whatever the list says. A run that records
spacing_source (a single-plane run, or io.voxel_size_um) is exported on the spacing it ran with, from its
image_summary, not on the file's.

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
# single-plane images (P2b): planar_helpers' 256 x 256 plane, and a plane of 1030 x 1100 px that is tiled
PLANE, PLANE_ID = "plane_001.nd2", "plane_001"
BIG, BIG_ID = "plane_big_001.nd2", "plane_big_001"
BIG_SHAPE = (1030, 1100)                 # tiles of 515 x 550 px: rows 0 and 515, columns 0 and 550
# disk centres of the big plane: none in its bottom-right tile, and one across the row-515 border on the left
BIG_DISKS = [(150, 150), (150, 400), (300, 700), (150, 950), (700, 200), (850, 400), (515, 275)]


def _plane_meta(data) -> dict:
    """An .nd2 header of a single plane: no Z axis (planar_helpers.nd2_meta, two_d)."""
    import planar_helpers as PH

    return PH.nd2_meta(data, (PH.DZ, PH.DY, PH.DX))


def _big_plane() -> np.ndarray:
    """(2, 1, 1030, 1100) uint16: bright disks of radius 8 px in the DAPI channel at BIG_DISKS."""
    rng = np.random.default_rng(5)
    yy, xx = np.indices(BIG_SHAPE)
    dna = np.full(BIG_SHAPE, 200.0)
    for cy, cx in BIG_DISKS:
        dna[(yy - cy) ** 2 + (xx - cx) ** 2 <= 8 ** 2] = 1400.0
    dna += rng.normal(0.0, 8.0, BIG_SHAPE)
    foci = 150.0 + rng.normal(0.0, 12.0, BIG_SHAPE)
    return np.clip(np.round(np.stack([dna, foci])), 0, 65535).astype(np.uint16)[:, None]


def _sources() -> dict:
    """{file name: (header, stack)} of the synthetic images, made once."""
    if not _STACKS:
        import planar_helpers as PH

        _STACKS[GONAD] = (GONAD_META, _synthetic_stack())
        stack, _planted = yeast_like_stack()
        c, z, y, x = stack.data.shape
        _STACKS[YEAST] = ({"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": stack.spacing, "spacing_ok": True,
                           "channel_names": list(stack.channel_names), "dtype": "float32", "is_2d": False}, stack)
        for name, data in ((PLANE, PH.plane_cyx()[:, None]), (BIG, _big_plane())):
            _STACKS[name] = (_plane_meta(data), Stack(data=data, spacing=(PH.DZ, PH.DY, PH.DX),
                                                      channel_names=list(PH.CHANNELS), path=Path(name),
                                                      spacing_ok=True))
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


def test_the_command_hands_on_a_list_as_a_tuple_and_auto_as_none(monkeypatch, tmp_path):
    """The plugin's in-process export calls export_training_pairs(orientations=<tuple>) and compares that call
    with the one `germquant export-training --orientations LIST` makes (napari-germquant
    tests/test_training.py::test_settings_reproduce_the_command_line_arguments), so a list goes on as the tuple
    of its items, as before; auto goes on as None."""
    seen = []

    def record(results_dir, out_dir, **kw):
        seen.append(kw["orientations"])
        return pd.DataFrame(columns=TR.MANIFEST_COLS)

    monkeypatch.setattr(TR, "export_training_pairs", record)
    for value in (None, "auto", "AUTO", "xy,yz", " xy , xz ,yz "):
        argv = ["export-training", str(tmp_path), "--out", str(tmp_path / "p")]
        assert cli.main(argv + ([] if value is None else ["--orientations", value])) == 0
    assert seen == [None, None, None, ("xy", "yz"), ("xy", "xz", "yz")]


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
    # the recorded plane is the one exported (the record, not a new pick, says what the segmenter saw)
    import shutil

    moved = tmp_path / "recorded"
    shutil.copytree(yeast, moved)
    for name in ("run_manifest.json", f"{YEAST_ID}__stages.json"):
        rec_file = moved / name
        data = json.loads(rec_file.read_text(encoding="utf-8"))
        data["segment_projection"]["z"] = Z_FOCUS + 2
        rec_file.write_text(json.dumps(data), encoding="utf-8")
    m = TR.export_training_pairs(moved, tmp_path / "recorded_pairs")
    assert int(m["slice_index"].iloc[0]) == Z_FOCUS + 2
    assert np.array_equal(tifffile.imread(tmp_path / "recorded_pairs" / m["image"].iloc[0]),
                          _u16(stack.data[0][Z_FOCUS + 2]))
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


# ---- a batch tree: the run folder has no run_manifest.json ---------------------------------------------------

def test_a_projected_batch_run_exports_the_channel_it_segmented(tmp_path, monkeypatch):
    """A batch writes one run_manifest.json, at its root; the run folder holds <id>__stages.json, which the
    geometry is read from. The proj pair is still the channel the segmenter saw: the recipe's transmitted role,
    here channel 1 of the file, from the config of the batch's manifest. A manifest above the folder that
    records another config_hash is not the run's, and the role falls back to dna as before."""
    from test_projected_run import SP

    from germquant import provenance

    stack, _ = yeast_like_stack()
    swapped = Stack(data=stack.data[::-1].copy(), spacing=SP, channel_names=["Cy3", "Bright"], path="x.nd2",
                    spacing_ok=True)
    c, z, y, x = swapped.data.shape
    meta = {"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": SP, "spacing_ok": True,
            "channel_names": ["Cy3", "Bright"], "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda path, **kw: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda path, **kw: swapped)
    cfg = recipe(**{"spots.enabled": False})
    root = tmp_path / "batch"
    prov = provenance.write_manifest(root, config_hash=cfg.hash, config=cfg.as_dict(), extra={})
    pipeline.process_image(YEAST, cfg, root / YEAST_ID, prov=prov)
    assert not (root / YEAST_ID / "run_manifest.json").exists()
    man = TR.export_training_pairs(root, tmp_path / "pairs")
    r = man.iloc[0].to_dict()
    assert (r["orientation"], int(r["slice_index"]), r["channel_role"], int(r["channel_index"])) == \
        ("proj", Z_FOCUS, "transmitted", 1)
    assert np.array_equal(tifffile.imread(tmp_path / "pairs" / r["image"]), _u16(swapped.data[1][Z_FOCUS]))
    monkeypatch.chdir(root / YEAST_ID)                 # an export from inside the run folder: "." has no parents
    r = TR.export_training_pairs(".", tmp_path / "here").iloc[0].to_dict()
    assert (r["channel_role"], int(r["channel_index"])) == ("transmitted", 1)
    other = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    other["config_hash"] = "0" * len(cfg.hash)
    (root / "run_manifest.json").write_text(json.dumps(other), encoding="utf-8")
    r = TR.export_training_pairs(root, tmp_path / "other").iloc[0].to_dict()
    assert (r["orientation"], r["channel_role"], int(r["channel_index"])) == ("proj", "dna", 0)


# ---- single-plane runs (P2b): the plane, tiled above 1024 px -------------------------------------------------

SWAP, SWAP_ID = "plane_swap_001.nd2", "plane_swap_001"


@pytest.fixture(scope="module")
def planar_results(tmp_path_factory, results):
    """The single-plane runs (config/profiles/generic_2d.yaml, classical, spots off): plane (256 x 256, its file's
    0.8 um), strided (xy stride 2), override (io.voxel_size_um 0.4) and big (1030 x 1100, tiled); and a 3D run of
    the synthetic gonad with io.voxel_size_um [0.6, 0.25, 0.25] (gonad_override)."""
    import planar_helpers as PH

    root = tmp_path_factory.mktemp("planar")
    with pytest.MonkeyPatch.context() as mp:
        _patch_readers(mp)
        off = {"spots.enabled": False}
        pipeline.process_image(PLANE, PH.generic_2d(**off), root / "plane", prov=None)
        pipeline.process_image(PLANE, PH.generic_2d(**off), root / "strided", prov=None, xy_stride=2)
        pipeline.process_image(PLANE, PH.generic_2d(**off, **{"io.voxel_size_um": 0.4}), root / "override", prov=None)
        pipeline.process_image(BIG, PH.generic_2d(**off), root / "big", prov=None)
        cfg = load_config(ROOT / "config" / "config.yaml")
        cfg._data["segmentation"]["nuclei"]["method"] = "classical"
        cfg.set("io.voxel_size_um", [0.6, 0.25, 0.25])
        pipeline.process_image(GONAD, cfg, root / "gonad_override", prov=None)
    return root


def test_plane_tiles():
    assert TR.TILE_PX == 1024
    assert TR.plane_tiles(256, 256) == [(0, 0, 256, 256)]
    assert TR.plane_tiles(1024, 1024) == [(0, 0, 1024, 1024)]           # a plane that fits is one tile
    assert TR.plane_tiles(1025, 10) == [(0, 0, 512, 10), (512, 0, 513, 10)]
    assert TR.plane_tiles(*BIG_SHAPE) == [(0, 0, 515, 550), (0, 550, 515, 550), (515, 0, 515, 550),
                                          (515, 550, 515, 550)]
    for h, w in ((1, 1), (2049, 3000), (5000, 1023), (4097, 4096), (3071, 7)):
        tiles = TR.plane_tiles(h, w)
        cover = np.zeros((h, w), np.int32)
        for y0, x0, th, tw in tiles:
            assert 0 < th <= 1024 and 0 < tw <= 1024
            cover[y0:y0 + th, x0:x0 + tw] += 1
        assert (cover == 1).all(), (h, w)                                # every pixel in exactly one tile
        ny, nx = -(-h // 1024), -(-w // 1024)
        assert len(tiles) == ny * nx                                     # the fewest tiles of at most 1024 px
        assert max(t[2] for t in tiles) - min(t[2] for t in tiles) <= 1  # near-equal sizes on each axis
        assert max(t[3] for t in tiles) - min(t[3] for t in tiles) <= 1
        assert tiles == sorted(tiles)                                    # row by row


def test_planar_run_exports_xy_tiles(tmp_path, planar_results, readers):
    import planar_helpers as PH

    plane = planar_results / "plane"
    rec = json.loads((plane / "run_manifest.json").read_text(encoding="utf-8"))
    assert rec["planar"] is True and rec["spacing_source"] == "file" and "segment_projection" not in rec
    labels = _labels(plane, PLANE_ID)
    assert labels.shape == (256, 256)                                   # a single-plane run writes YX labels
    man = TR.export_training_pairs(plane, tmp_path / "plane")
    assert list(man.columns) == TR.MANIFEST_COLS + ["tile_y0", "tile_x0"]
    assert len(man) == 1
    r = man.iloc[0].to_dict()
    stem = f"{PLANE_ID}_xy_0000.tif"
    assert (r["orientation"], int(r["slice_index"]), r["image"], r["labels"]) == ("xy", 0, f"image/{stem}",
                                                                                f"labels/{stem}")
    assert (int(r["tile_y0"]), int(r["tile_x0"]), int(r["height"]), int(r["width"])) == (0, 0, 256, 256)
    assert (r["channel_role"], int(r["channel_index"]), float(r["pixel_um"]), float(r["z_zoom"])) == \
        ("dna", 0, 0.8, 1.0)
    assert int(r["n_labels"]) == len(np.unique(labels)) - 1 == len(PH.planted()["nuclei"])
    assert np.array_equal(tifffile.imread(tmp_path / "plane" / r["image"]), _u16(PH.plane_cyx()[0]))
    assert np.array_equal(tifffile.imread(tmp_path / "plane" / r["labels"]), labels)
    exp = _export_json(tmp_path / "plane")
    assert exp["orientations"] == ["xy"]
    assert exp["run_orientations"] == [{"image_id": PLANE_ID, "geometry": "planar", "orientations": ["xy"],
                                        "tile_px": 1024, "spacing_source": "file"}]
    # an explicit list names the slices of 3D runs: a single-plane run still gives its plane
    listed = TR.export_training_pairs(plane, tmp_path / "listed", orientations=("xz",))
    pd.testing.assert_frame_equal(listed, man)
    assert _export_json(tmp_path / "listed")["orientations"] == ["xy"]
    # read again at the run's geometry: xy stride 2 halves the plane and doubles the pixel size
    s = TR.export_training_pairs(planar_results / "strided", tmp_path / "strided").iloc[0]
    assert (int(s["height"]), int(s["width"]), float(s["pixel_um"])) == (128, 128, 1.6)
    assert np.array_equal(tifffile.imread(tmp_path / "strided" / s["image"]), _u16(PH.plane_cyx()[0][::2, ::2]))
    # io.voxel_size_um: the pixel size the run used (its image_summary), not the file's 0.8 um
    o = TR.export_training_pairs(planar_results / "override", tmp_path / "override").iloc[0]
    assert float(o["pixel_um"]) == 0.4
    assert _export_json(tmp_path / "override")["run_orientations"][0]["spacing_source"] == "override"
    # the big plane: the fewest tiles of at most 1024 px; the empty bottom-right tile (number 3) is not written
    big = planar_results / "big"
    big_labels = _labels(big, BIG_ID)
    assert big_labels.shape == BIG_SHAPE
    b = TR.export_training_pairs(big, tmp_path / "big")
    assert [int(v) for v in b["slice_index"]] == [0, 1, 2]
    assert [(int(y), int(x)) for y, x in zip(b["tile_y0"], b["tile_x0"])] == [(0, 0), (0, 550), (515, 0)]
    assert set(zip(b["height"].astype(int), b["width"].astype(int))) == {(515, 550)}
    data = _sources()[BIG][1].data
    seen: dict[int, set] = {}
    for _, row in b.iterrows():
        y0, x0, h, w = int(row["tile_y0"]), int(row["tile_x0"]), int(row["height"]), int(row["width"])
        im = tifffile.imread(tmp_path / "big" / row["image"])
        lab = tifffile.imread(tmp_path / "big" / row["labels"])
        assert np.array_equal(im, _u16(data[0, 0, y0:y0 + h, x0:x0 + w]))
        assert np.array_equal(lab, big_labels[y0:y0 + h, x0:x0 + w]) and lab.dtype == np.int32
        ids = set(np.unique(lab)) - {0}
        assert int(row["n_labels"]) == len(ids)
        assert row["image"] == f"image/{BIG_ID}_xy_{int(row['slice_index']):04d}.tif"
        for i in ids:
            seen.setdefault(int(i), set()).add(int(row["slice_index"]))
    assert not big_labels[515:, 550:].any()                              # the tile that was left out is empty
    assert set(seen) == set(np.unique(big_labels)) - {0}                 # every object is in some tile
    across = int(big_labels[515, 275])
    assert across and seen[across] == {0, 2}                             # the disk across the border: both pieces
    assert _export_json(tmp_path / "big")["run_orientations"][0]["tile_px"] == 1024
    assert cli.main(["export-training", str(big), "--out", str(tmp_path / "big_cli")]) == 0
    assert _manifest_text(tmp_path / "big_cli") == _manifest_text(tmp_path / "big")
    # the pairs load with their tile number as the slice index; the card calls them 2D, also without manifest.csv
    pairs = TR.load_pairs(tmp_path / "big")
    assert [(p["orientation"], p["slice_index"]) for p in pairs] == [("xy", 0), ("xy", 1), ("xy", 2)]
    assert TR._card_fields(tmp_path / "big", pairs)["dims"] == "2d"
    (tmp_path / "big" / "manifest.csv").unlink()
    named = TR.load_pairs(tmp_path / "big")
    assert [(p["image_id"], p["slice_index"]) for p in named] == [(BIG_ID, 0), (BIG_ID, 1), (BIG_ID, 2)]
    assert TR._card_fields(tmp_path / "big", named)["dims"] == "2d"


def test_a_tree_of_3d_and_single_plane_runs(tmp_path, results, planar_results, readers):
    import shutil

    tree = tmp_path / "tree"
    shutil.copytree(results / "mixed" / "gonad", tree / "gonad")
    shutil.copytree(planar_results / "plane", tree / "plane")
    man = TR.export_training_pairs(tree, tmp_path / "pairs")
    exp = _export_json(tmp_path / "pairs")
    assert exp["orientations"] == ["xy", "xz", "yz"]                      # xy is in the slice list already
    assert [(e["image_id"], e["geometry"], e["orientations"]) for e in exp["run_orientations"]] == [
        (GONAD_ID, "3d", ["xy", "xz", "yz"]), (PLANE_ID, "planar", ["xy"])]
    gonad_rows = man[man["image_id"] == GONAD_ID]
    assert gonad_rows["tile_y0"].isna().all() and len(gonad_rows) == 9
    # the 3D run's rows are the ones it gives alone, with the two tile cells empty
    TR.export_training_pairs(tree / "gonad", tmp_path / "alone")
    lines = _manifest_text(tmp_path / "pairs").splitlines()
    assert lines[0] == HEADER + ",tile_y0,tile_x0"
    alone = _manifest_text(tmp_path / "alone").splitlines()[1:]
    assert [ln for ln in lines if ln.startswith(GONAD_ID + ",")] == [ln + ",," for ln in alone]
    plane_rows = man[man["image_id"] == PLANE_ID]
    assert list(plane_rows["orientation"]) == ["xy"] and int(plane_rows["tile_y0"].iloc[0]) == 0
    assert TR._card_fields(tmp_path / "pairs", TR.load_pairs(tmp_path / "pairs"))["dims"] == "2d+3d"
    # an explicit list: the 3D slices it names, then the single plane's xy
    TR.export_training_pairs(tree, tmp_path / "yz", orientations=("yz",))
    assert _export_json(tmp_path / "yz")["orientations"] == ["yz", "xy"]


def test_a_run_with_the_voxel_size_override_is_exported_on_its_spacing(tmp_path, results, planar_results, readers):
    """io.voxel_size_um (D6) is applied after the read, so the export takes the run's spacing from its
    image_summary: the xz and yz slices are stretched by the run's dz / dy and pixel_um is the run's dy."""
    run = planar_results / "gonad_override"
    rec = json.loads((run / "run_manifest.json").read_text(encoding="utf-8"))
    assert rec["spacing_source"] == "override"
    man = TR.export_training_pairs(run, tmp_path / "pairs")
    assert list(man.columns) == TR.MANIFEST_COLS                       # a 3D run: no tile columns
    assert set(man["pixel_um"].astype(float)) == {0.25}
    side = man[man["orientation"] != "xy"]
    assert set(side["z_zoom"].astype(float)) == {0.6 / 0.25} and set(side["height"].astype(int)) == {48}
    stack, labels = _synthetic_stack(), _labels(run, GONAD_ID)
    for _, r in side.iterrows():
        im2, lab2, px, z_zoom = TR.isotropic_slice(stack.data[0], labels, r["orientation"], int(r["slice_index"]),
                                                   (0.6, 0.25, 0.25))
        assert np.array_equal(tifffile.imread(tmp_path / "pairs" / r["image"]), _u16(im2))
        assert np.array_equal(tifffile.imread(tmp_path / "pairs" / r["labels"]), lab2)
    # the same image without the override: the file's spacing, as before
    base = TR.export_training_pairs(results / "mixed" / "gonad", tmp_path / "base")
    assert set(base[base["orientation"] != "xy"]["z_zoom"].astype(float)) == {2.0}
    assert set(base["pixel_um"].astype(float)) == {0.2}


def test_a_planar_batch_run_exports_the_channel_it_segmented(tmp_path, monkeypatch):
    """A batch run folder holds no run_manifest.json; its single-plane record sends the export to the batch's
    manifest for the segmented role (here foci, channel 1 of the file, which holds the nuclei), as for a
    projected run of a batch (review S1 of part 1)."""
    import planar_helpers as PH

    from germquant import provenance

    _patch_readers(monkeypatch)
    swapped = PH.plane_cyx()[::-1].copy()[:, None]           # channel 0 (DAPI) flat, channel 1 (foci) the nuclei
    _STACKS[SWAP] = (_plane_meta(swapped), Stack(data=swapped, spacing=(PH.DZ, PH.DY, PH.DX),
                                                 channel_names=list(PH.CHANNELS), path=Path(SWAP), spacing_ok=True))
    cfg = PH.generic_2d(**{"spots.enabled": False, "segmentation.nuclei.channel": "foci"})
    root = tmp_path / "batch"
    prov = provenance.write_manifest(root, config_hash=cfg.hash, config=cfg.as_dict(), extra={})
    pipeline.process_image(SWAP, cfg, root / SWAP_ID, prov=prov)
    assert not (root / SWAP_ID / "run_manifest.json").exists()
    r = TR.export_training_pairs(root, tmp_path / "pairs").iloc[0].to_dict()
    assert (r["orientation"], r["channel_role"], int(r["channel_index"])) == ("xy", "foci", 1)
    assert np.array_equal(tifffile.imread(tmp_path / "pairs" / r["image"]), _u16(swapped[1, 0]))
