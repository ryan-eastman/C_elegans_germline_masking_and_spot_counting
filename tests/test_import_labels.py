"""import-labels: labels made elsewhere, written as training pairs (docs/TRAINING.md, step 3; docs/DESIGN_system_packs.md
sections 6.6 and 10 step 3; P5 acceptance 9; decisions D23 and D28 of the Train B plan).

``germquant import-labels SRC --format tiff|png|cellacdc-npz|cellpose-npy --image IMG --pairs DIR`` reads a label
image made in another program and writes the export layout: image/ and labels/ TIFFs and manifest.csv rows whose
source is ``import:<tool>``. Labels of one plane give the plane (tiled above 1024 px, as a single-plane run's export),
labels of every plane give xy, xz and yz slices drawn as export-training draws them. A Cellpose ``_seg.npy`` is a
pickle and is read only with --allow-pickle, with a warning; an image that a committed benchmark split
(benchmarks/ext/<dataset>/split.json with converted.json) holds as test data is refused by its file sha256 or its
pixels_sha256; with no committed split at all the import warns "no committed splits found". In a fresh
interpreter import_labels imports germquant.pipeline (torch) before it reads the image (D28).

The images reach the engine through the .nd2 reader seam of process_image (pipeline.read_nd2_metadata and
read_stack, as tests/test_training_export_geometry.py feeds its runs), except in the fresh-interpreter test, which
reads an OME-TIFF through bioio.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile
from packs_helpers import subprocess_env

from germquant import cli, config, pipeline
from germquant.io import labels as LABELS
from germquant.io.nd2_reader import Stack
from germquant.io.pixels import pixels_sha256
from germquant.segment import training as TR

ROOT = Path(__file__).resolve().parents[1]
CHANNELS = ["DAPI", "foci"]
SP2 = (1.0, 0.65, 0.65)          # dz (no z step: a one-plane file), dy, dx um
SP3 = (0.5, 0.2, 0.2)
SOURCE_IDS = (5, 9, 12)          # the ids of the labels files, relabelled 1, 2, 3 on import


def _disks(shape, centres, r, ids) -> np.ndarray:
    yy, xx = np.indices(shape[-2:])
    lab = np.zeros(shape[-2:], np.uint16)
    for (cy, cx), i in zip(centres, ids):
        lab[(yy - cy) ** 2 + (xx - cx) ** 2 <= r * r] = i
    return lab


def _labels_2d() -> np.ndarray:
    """(96, 120) uint16 labels with the ids 5, 9 and 12."""
    return _disks((96, 120), [(30, 30), (40, 80), (70, 50)], 11, SOURCE_IDS)


def _labels_3d() -> np.ndarray:
    """(6, 48, 64) uint16: three blobs through planes 1 to 4, with the ids 5, 9 and 12."""
    lab = np.zeros((6, 48, 64), np.uint16)
    plane = _disks((48, 64), [(12, 14), (24, 40), (36, 20)], 7, SOURCE_IDS)
    lab[1:5] = plane
    return lab


def _pixels(lab: np.ndarray, seed: int) -> np.ndarray:
    """Two uint16 channels: bright objects on a dim field in channel 0, a flat noisy channel 1."""
    rng = np.random.default_rng(seed)
    dna = 200.0 + 1200.0 * (lab > 0) + rng.normal(0, 8, lab.shape)
    foci = 150.0 + rng.normal(0, 12, lab.shape)
    return np.clip(np.round(np.stack([dna, foci])), 0, 65535).astype(np.uint16)


class Images:
    """The images of these tests, behind the .nd2 reader seam: {name: (header, data per timepoint)}."""

    def __init__(self, folder: Path):
        self.folder = folder
        self.src: dict[str, tuple[dict, list[np.ndarray], tuple]] = {}
        self.reads: list[tuple] = []

    def add(self, name: str, data, spacing, *, two_d: bool, series: int = 1) -> Path:
        frames = [np.asarray(d) for d in (data if isinstance(data, list) else [data])]
        c, z, y, x = frames[0].shape
        sizes = ({"T": series} if series > 1 else {}) | {"C": c} | ({} if two_d else {"Z": z}) | {"Y": y, "X": x}
        meta = {"sizes": sizes, "spacing": tuple(spacing), "spacing_ok": True, "channel_names": list(CHANNELS),
                "dtype": "uint16", "is_2d": two_d, "n_timepoints": series, "period_s": None, "t_seconds": None}
        self.src[name] = (meta, frames, tuple(spacing))
        path = self.folder / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"stand-in " + name.encode("ascii") + b" " + frames[0].tobytes()[:64])
        return path

    def data(self, name: str, t: int = 0) -> np.ndarray:
        return self.src[name][1][t]

    def patch(self, mp) -> None:
        def read_meta(path, t=None):
            return dict(self.src[Path(str(path)).name][0])

        def read_stack(path, *, xy_stride=1, z_range=None, as_float=False, t=None):
            name = Path(str(path)).name
            self.reads.append((name, xy_stride, z_range, t))
            meta, frames, sp = self.src[name]
            arr = frames[0 if t is None else int(t)]
            return Stack(data=np.ascontiguousarray(arr), spacing=sp, channel_names=list(CHANNELS), path=Path(path),
                         spacing_ok=True, t_index=t)

        mp.setattr(pipeline, "read_nd2_metadata", read_meta)
        mp.setattr(pipeline, "read_stack", read_stack)


@pytest.fixture
def images(tmp_path, monkeypatch):
    im = Images(tmp_path / "images")
    im.patch(monkeypatch)
    # no committed split anywhere: the config roots are an empty folder (the checkout's benchmarks/ext/ may come)
    empty = tmp_path / "no_splits_root"
    empty.mkdir()
    monkeypatch.setattr(config, "config_roots", lambda: [empty])
    return im


def _manifest(pairs: Path) -> pd.DataFrame:
    return pd.read_csv(pairs / "manifest.csv", low_memory=False)


def _write_labels(path: Path, kind: str, lab: np.ndarray) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if kind == "tiff":
        tifffile.imwrite(path, lab)
    elif kind == "png":
        from PIL import Image

        Image.fromarray(lab.astype(np.uint16)).save(path)
    elif kind == "npz":
        np.savez_compressed(path, lab)
    elif kind == "npy":
        np.save(path, np.array({"masks": lab, "outlines": np.zeros_like(lab)}, dtype=object), allow_pickle=True)
    return path


def _relabelled(lab: np.ndarray) -> np.ndarray:
    out = np.zeros(lab.shape, np.int32)
    for new, old in enumerate(SOURCE_IDS, start=1):
        out[lab == old] = new
    return out


def _u16(a) -> np.ndarray:
    return np.clip(np.rint(np.asarray(a, dtype=np.float64)), 0, 65535).astype(np.uint16)


# ---- the formats (acceptance 9) -----------------------------------------------------------------------------

CASES = {
    # case: (format, file name, the labels as written, image kind, frame, tool)
    "tiff": ("tiff", "cells_001_masks.tif", "2d", 0, None),
    "png": ("png", "cells_001_masks.png", "2d", 0, None),
    "npz2d": ("cellacdc-npz", "cells_001_segm.npz", "2d", 0, None),
    "npz3d": ("cellacdc-npz", "stack_001_segm.npz", "3d", 0, None),
    "npzT": ("cellacdc-npz", "stack_001_segmT.npz", "3d", 1, "cellacdc"),
}


def _case_files(images: Images, tmp_path: Path, case: str):
    fmt, fname, dims, frame, _tool = CASES[case]
    if dims == "2d":
        lab = _labels_2d()
        img = images.add("cells_001.nd2", _pixels(lab, 1)[:, None], SP2, two_d=True)
    else:
        lab = _labels_3d()
        img = images.add("stack_001.nd2", _pixels(lab, 2), SP3, two_d=False)
    kind = {"tiff": "tiff", "png": "png", "cellacdc-npz": "npz"}[fmt]
    written = lab
    if case == "npzT":
        # (T, Z, Y, X): frame 1 holds the labels, frame 0 something else
        other = np.zeros_like(lab)
        other[2:4, 30:40, 40:60] = 3
        written = np.stack([other, lab])
    src = _write_labels(tmp_path / "labels" / fname, kind, written)
    return fmt, src, img, lab, frame


@pytest.mark.parametrize("case", list(CASES))
def test_formats(images, tmp_path, case):
    fmt, src, img, lab, frame = _case_files(images, tmp_path, case)
    from germquant.import_labels import import_labels

    pairs = tmp_path / "pairs"
    res = import_labels(src, fmt, img, pairs, frame=frame)
    tool = {"tiff": "tiff", "png": "png", "cellacdc-npz": "cellacdc"}[fmt]
    stem = Path(img).name[:-len(".nd2")]
    data = images.data(Path(img).name)
    assert res["image_id"] == stem and res["series_id"] == stem and res["t_index"] is None
    assert res["source"] == f"import:{tool}" and res["tool"] == tool and res["format"] == fmt
    assert res["n_objects"] == 3 and res["channel_index"] == 0 and res["channel_name"] == "DAPI"
    assert res["pixels_sha256"] == pixels_sha256(data)
    assert res["image_sha256"] == hashlib.sha256(Path(img).read_bytes()).hexdigest()
    assert res["labels_sha256"] == hashlib.sha256(Path(src).read_bytes()).hexdigest()
    assert any("no committed splits found" in w for w in res["warnings"])
    assert any("channel 0" in w and "DAPI" in w for w in res["warnings"])     # two channels, none named
    man = _manifest(pairs)
    two_d = CASES[case][2] == "2d"
    # the import columns after MANIFEST_COLS: the source, and the tile offsets of a plane (2D labels only)
    assert list(man.columns) == TR.MANIFEST_COLS + ["source"] + (["tile_y0", "tile_x0"] if two_d else [])
    assert set(man["source"]) == {f"import:{tool}"} and res["n_pairs"] == len(man) == len(res["pairs"])
    assert set(man["image_id"]) == {stem} and set(man["channel_index"]) == {0}
    assert man["channel_role"].isna().all()                        # the role the labels were drawn for is unknown
    assert man["source_file"].iloc[0] == os.path.abspath(img)
    want = _relabelled(lab)
    if two_d:
        assert res["dims"] == "2d" and len(man) == 1
        r = man.iloc[0]
        assert (r["orientation"], int(r["slice_index"]), int(r["tile_y0"]), int(r["tile_x0"])) == ("xy", 0, 0, 0)
        assert (int(r["height"]), int(r["width"]), float(r["pixel_um"]), float(r["z_zoom"])) == (96, 120, 0.65, 1.0)
        assert int(r["n_labels"]) == 3 and r["image"] == f"image/{stem}_xy_0000.tif"
        assert np.array_equal(tifffile.imread(pairs / r["image"]), _u16(data[0, 0]))
        got = tifffile.imread(pairs / r["labels"])
        assert got.dtype == np.int32 and np.array_equal(got, want)           # 1..N in the order of the file's ids
    else:
        assert res["dims"] == "3d"
        # the slices export-training draws for a 3D run: seeded per image, 3 per orientation, xz and yz isotropic
        rng = TR._rng_for(0, stem)
        expect = []
        for orient in TR.ORIENTATIONS:
            for i in TR.choose_slices(TR.slice_counts(want, orient), 3, 3, rng):
                expect.append((orient, i))
        assert [(r["orientation"], int(r["slice_index"])) for _, r in man.iterrows()] == expect
        for _, r in man.iterrows():
            im2, lab2, px, z_zoom = TR.isotropic_slice(data[0], want, r["orientation"], int(r["slice_index"]), SP3)
            assert np.array_equal(tifffile.imread(pairs / r["image"]), _u16(im2))
            assert np.array_equal(tifffile.imread(pairs / r["labels"]), lab2)
            assert (float(r["pixel_um"]), float(r["z_zoom"])) == (0.2, 1.0 if r["orientation"] == "xy" else 2.5)
    # the pairs load for training, with their orientation, under the image's id
    loaded = TR.load_pairs(pairs)
    assert [(p["image_id"], p["orientation"]) for p in loaded] == [(stem, o) for o in man["orientation"]]
    # read at full resolution, every plane, no timepoint for a fixed image
    assert images.reads[-1] == (Path(img).name, 1, None, None)


def test_a_series_image_is_read_at_the_frame_and_named_by_its_timepoint(images, tmp_path):
    """A (T, Y, X) segm.npz of a one-plane time series: --frame picks the labels' frame and the image's
    timepoint, and the pairs carry the timepoint's image id (<stem>_t<NNN>) in the series' group. A labels file
    without a T axis, or with a T axis of one frame, holds the labels of the timepoint --frame names; a warning
    names the timepoint the pairs take."""
    from germquant.import_labels import import_labels

    lab = _labels_2d()
    frames = [_pixels(np.roll(lab, 3 * t, axis=1), 10 + t)[:, None] for t in range(3)]
    img = images.add("movie_001.nd2", frames, SP2, two_d=True, series=3)
    stackT = np.stack([np.roll(lab, 3 * t, axis=1) for t in range(3)])
    src = _write_labels(tmp_path / "movie_001_segm.npz", "npz", stackT)
    res = import_labels(src, "cellacdc-npz", img, tmp_path / "pairs", frame=2)
    assert (res["image_id"], res["series_id"], res["t_index"]) == ("movie_001_t002", "movie_001", 2)
    assert images.reads[-1] == ("movie_001.nd2", 1, None, 2)
    r = _manifest(tmp_path / "pairs").iloc[0]
    assert (r["image_id"], r["series_id"], int(r["t_index"])) == ("movie_001_t002", "movie_001", 2)
    assert np.array_equal(tifffile.imread(tmp_path / "pairs" / r["image"]), _u16(frames[2][0, 0]))
    assert np.array_equal(tifffile.imread(tmp_path / "pairs" / r["labels"]), _relabelled(np.roll(lab, 6, axis=1)))
    assert res["pixels_sha256"] == pixels_sha256(frames[2])
    # the timepoint the pairs take is named, also when --frame is left at 0
    assert any("time series of 3 timepoints" in w and "timepoint 2" in w for w in res["warnings"])
    # labels of one timepoint saved without a T axis: --frame names the timepoint they belong to
    one = _write_labels(tmp_path / "movie_001_t1_masks.tif", "tiff", np.roll(lab, 3, axis=1))
    res1 = import_labels(one, "tiff", img, tmp_path / "pairs_one", frame=1)
    assert (res1["image_id"], res1["t_index"]) == ("movie_001_t001", 1)
    r1 = _manifest(tmp_path / "pairs_one").iloc[0]
    assert np.array_equal(tifffile.imread(tmp_path / "pairs_one" / r1["image"]), _u16(frames[1][0, 0]))
    assert np.array_equal(tifffile.imread(tmp_path / "pairs_one" / r1["labels"]),
                          _relabelled(np.roll(lab, 3, axis=1)))
    # left at frame 0, the same file is timepoint 0's, and the warning names the timepoint the pairs took
    res0 = import_labels(one, "tiff", img, tmp_path / "pairs_zero")
    assert res0["t_index"] == 0 and any("timepoint 0" in w for w in res0["warnings"])
    # a T axis of one frame (a (1, 1, Y, X) segm.npz) is one timepoint's labels too
    single = _write_labels(tmp_path / "movie_001_t1_segm.npz", "npz", np.roll(lab, 3, axis=1)[None, None])
    res_single = import_labels(single, "cellacdc-npz", img, tmp_path / "pairs_single", frame=1)
    assert (res_single["image_id"], res_single["t_index"]) == ("movie_001_t001", 1)
    r_single = _manifest(tmp_path / "pairs_single").iloc[0]
    assert np.array_equal(tifffile.imread(tmp_path / "pairs_single" / r_single["labels"]),
                          _relabelled(np.roll(lab, 3, axis=1)))
    # a frame the image does not have is refused before anything is written
    with pytest.raises(ValueError):
        import_labels(src, "cellacdc-npz", img, tmp_path / "pairs_bad", frame=3)
    with pytest.raises(ValueError):
        import_labels(one, "tiff", img, tmp_path / "pairs_bad", frame=3)
    assert not (tmp_path / "pairs_bad").exists()


def test_seg_npy_needs_allow_pickle(images, tmp_path, caplog, capsys):
    from germquant.import_labels import import_labels

    lab = _labels_2d()
    img = images.add("cells_001.nd2", _pixels(lab, 1)[:, None], SP2, two_d=True)
    src = _write_labels(tmp_path / "cells_001_seg.npy", "npy", lab)
    pairs = tmp_path / "pairs"
    with pytest.raises(LABELS.LabelsError, match="pickle"):
        import_labels(src, "cellpose-npy", img, pairs)
    assert images.reads == [] and not pairs.exists()          # refused before the image or the file was read
    argv = ["import-labels", str(src), "--format", "cellpose-npy", "--image", str(img), "--pairs", str(pairs)]
    assert cli.main(argv) == 2
    err = capsys.readouterr().err
    assert "REFUSED:" in err and "--allow-pickle" in err and not pairs.exists()
    assert cli.main(argv + ["--json"]) == 2
    out = json.loads(capsys.readouterr().out)
    assert out["kind"] == "LabelsError" and "pickle" in out["error"]
    # with --allow-pickle the masks are read, with a warning
    caplog.clear()
    with caplog.at_level("WARNING", logger="germquant.import_labels"):
        res = import_labels(src, "cellpose-npy", img, pairs, allow_pickle=True)
    assert res["source"] == "import:cellpose" and res["n_pairs"] == 1
    assert any("pickle" in w and "can run code" in w for w in res["warnings"])
    assert any("pickle" in r.getMessage() for r in caplog.records)
    got = tifffile.imread(pairs / _manifest(pairs)["labels"].iloc[0])
    assert np.array_equal(got, _relabelled(lab))


# ---- the test-split refusal (design 6.6, D23) ---------------------------------------------------------------

def _split_root(root: Path, dataset: str, items: dict, test=(), tune=(), val=()) -> Path:
    ext = root / "benchmarks" / "ext" / dataset
    ext.mkdir(parents=True, exist_ok=True)
    (ext / "split.json").write_text(json.dumps({"format": 1, "seed": 0, "tune": list(tune), "val": list(val),
                                                "test": list(test)}), encoding="ascii")
    (ext / "converted.json").write_text(json.dumps({"format": 1, "items": items}), encoding="ascii")
    return root


def test_test_split_image_refused(images, tmp_path, monkeypatch, capsys):
    from germquant.import_labels import SplitTestImageError, import_labels

    lab = _labels_2d()
    img = images.add("cells_001.nd2", _pixels(lab, 1)[:, None], SP2, two_d=True)
    src = _write_labels(tmp_path / "cells_001_masks.tif", "tiff", lab)
    file_sha = hashlib.sha256(img.read_bytes()).hexdigest()
    pix_sha = pixels_sha256(images.data("cells_001.nd2"))
    other = {"file": "x.ome.tif", "file_sha256": "1" * 64, "pixels_sha256": "2" * 64}
    pairs = tmp_path / "pairs"

    def roots(*paths):
        monkeypatch.setattr(config, "config_roots", lambda: [Path(p) for p in paths])

    # held out as test by its file sha256: refused, nothing written
    by_file = _split_root(tmp_path / "r1", "ds_one", {"t1": {**other, "file_sha256": file_sha}, "a": other},
                          test=["t1"], tune=["a"])
    roots(tmp_path / "no_splits_root", by_file)
    with pytest.raises(SplitTestImageError, match="ds_one") as e:
        import_labels(src, "tiff", img, pairs)
    assert "t1" in str(e.value) and "file sha256" in str(e.value)
    assert not pairs.exists()
    # by its pixels: the same pixels in another container are the same test image
    by_pixels = _split_root(tmp_path / "r2", "ds_two", {"t9": {**other, "pixels_sha256": pix_sha}}, test=["t9"])
    roots(by_pixels)
    with pytest.raises(SplitTestImageError, match="pixels_sha256"):
        import_labels(src, "tiff", img, pairs)
    argv = ["import-labels", str(src), "--format", "tiff", "--image", str(img), "--pairs", str(pairs), "--json"]
    assert cli.main(argv) == 2
    out = json.loads(capsys.readouterr().out)
    assert out["kind"] == "SplitTestImageError" and "ds_two" in out["error"] and not pairs.exists()
    # the same image in a tune or val partition is training material: imported, the splits named as checked
    tune_only = _split_root(tmp_path / "r3", "ds_three", {"a": {**other, "file_sha256": file_sha},
                                                         "b": {**other, "pixels_sha256": pix_sha}},
                            tune=["a"], val=["b"])
    roots(tune_only, by_file / "not_there")
    res = import_labels(src, "tiff", img, pairs)
    assert res["n_pairs"] == 1 and res["splits_checked"] == ["ds_three"]
    assert not any("no committed splits found" in w for w in res["warnings"])
    # a split whose test items carry no digests cannot be checked: said, not passed silently
    blind = _split_root(tmp_path / "r4", "ds_blind", {}, test=["z1", "z2"])
    roots(blind)
    res = import_labels(src, "tiff", img, tmp_path / "pairs_blind")
    assert any("ds_blind" in w and "z1" in w for w in res["warnings"])


def test_no_committed_splits_warns(images, tmp_path, caplog):
    from germquant.import_labels import import_labels

    lab = _labels_2d()
    img = images.add("cells_001.nd2", _pixels(lab, 1)[:, None], SP2, two_d=True)
    src = _write_labels(tmp_path / "m.tif", "tiff", lab)
    with caplog.at_level("WARNING", logger="germquant.import_labels"):
        res = import_labels(src, "tiff", img, tmp_path / "pairs", channel="DAPI")
    assert res["splits_checked"] == []
    assert any("no committed splits found" in w for w in res["warnings"])
    assert any("no committed splits found" in r.getMessage() for r in caplog.records)
    assert not any("channel 0" in w for w in res["warnings"])        # the channel was named


def test_pixels_sha256_is_stable():
    a = np.arange(12, dtype=np.uint16).reshape(3, 4)
    want = hashlib.sha256(b"<u2|3,4|" + a.astype("<u2").tobytes()).hexdigest()
    assert pixels_sha256(a) == want
    # pinned: P7's converted.json records this digest, so it must never change for the same pixels
    assert want == "4920aaab215d4d6c539b2f372b7068d52a85cd709277f2bd4edea5361a53d28f"
    assert pixels_sha256(a.astype(">u2")) == want                # the byte order of the array does not count
    assert pixels_sha256(np.asfortranarray(a)) == want            # nor its memory layout: C order is hashed
    big = np.arange(24, dtype=np.uint16).reshape(3, 8)
    assert pixels_sha256(big[:, ::2]) == pixels_sha256(np.ascontiguousarray(big[:, ::2]))
    assert pixels_sha256(a.astype(np.int32)) != want             # the dtype and the shape do
    assert pixels_sha256(a.reshape(4, 3)) != want
    assert pixels_sha256(a.reshape(1, 3, 4)) != want
    assert pixels_sha256(np.zeros((0, 4), np.uint8)) == hashlib.sha256(b"|u1|0,4|").hexdigest()
    with pytest.raises(ValueError):
        pixels_sha256(np.array([object()], dtype=object))


# ---- refusals ----------------------------------------------------------------------------------------------

def test_refusals_before_anything_is_written(images, tmp_path, capsys):
    from germquant.import_labels import ImportLabelsError, import_labels

    lab2, lab3 = _labels_2d(), _labels_3d()
    img2 = images.add("cells_001.nd2", _pixels(lab2, 1)[:, None], SP2, two_d=True)
    img3 = images.add("stack_001.nd2", _pixels(lab3, 2), SP3, two_d=False)
    tif2 = _write_labels(tmp_path / "l" / "m2.tif", "tiff", lab2)
    png2 = _write_labels(tmp_path / "l" / "m2.png", "png", lab2)
    npz3 = _write_labels(tmp_path / "l" / "m3.npz", "npz", lab3)
    flat3 = _write_labels(tmp_path / "l" / "m3_plane.tif", "tiff", lab3[2])   # one plane on the stack's own grid
    zyx = tmp_path / "l" / "m3.ome.tif"                       # a z stack by its axes, for an image of one plane
    tifffile.imwrite(zyx, lab3, ome=True, metadata={"axes": "ZYX"})
    pairs = tmp_path / "pairs"
    bad = [
        (dict(src=tif2, fmt="jpeg", image=img2), ValueError, "format"),
        (dict(src=tif2, fmt="png", image=img2), ImportLabelsError, "png"),                 # the name says TIFF
        (dict(src=npz3, fmt="tiff", image=img3), ImportLabelsError, "tiff"),
        (dict(src=flat3, fmt="tiff", image=img3), ImportLabelsError, r"holds one plane \(2D labels\)"),
        (dict(src=zyx, fmt="tiff", image=img2), LABELS.LabelsError, "planes"),              # 6 planes, image 1
        (dict(src=tif2, fmt="tiff", image=img2, frame=1), LABELS.LabelsError, "--frame 1 names a frame"),  # no T
        (dict(src=png2, fmt="png", image=img2, frame=1), LABELS.LabelsError, "--frame 1 names a frame"),
        (dict(src=tif2, fmt="tiff", image=img2, tool="my tool"), ValueError, "tool"),
        (dict(src=tif2, fmt="tiff", image=img2, channel=2), ImportLabelsError, "channel"),
        (dict(src=tif2, fmt="tiff", image=img2, channel="Cy5"), ImportLabelsError, "Cy5"),
        (dict(src=tmp_path / "l" / "gone.tif", fmt="tiff", image=img2), FileNotFoundError, "gone"),
        (dict(src=tif2, fmt="tiff", image=tmp_path / "images" / "gone.nd2"), FileNotFoundError, "gone"),
    ]
    for kw, exc, text in bad:
        with pytest.raises(exc, match=text):
            import_labels(kw.pop("src"), kw.pop("fmt"), kw.pop("image"), pairs, **kw)
        assert not pairs.exists(), (exc, text)
    # labels on another grid than the image; labels without any object
    small = _write_labels(tmp_path / "l" / "small.tif", "tiff", lab2[:90, :100])
    with pytest.raises(ImportLabelsError, match=r"\(90, 100\)"):
        import_labels(small, "tiff", img2, pairs)
    empty = _write_labels(tmp_path / "l" / "empty.tif", "tiff", np.zeros_like(lab2))
    with pytest.raises(ImportLabelsError, match="no object"):
        import_labels(empty, "tiff", img2, pairs)
    assert not pairs.exists()
    # the command: exit 2 and one REFUSED line; with --json the error and its kind
    argv = ["import-labels", str(flat3), "--format", "tiff", "--image", str(img3), "--pairs", str(pairs)]
    assert cli.main(argv) == 2
    refused = [ln for ln in capsys.readouterr().err.splitlines() if ln.startswith("REFUSED:")]
    assert len(refused) == 1 and "holds one plane (2D labels)" in refused[0]
    assert cli.main(argv + ["--json"]) == 2
    assert json.loads(capsys.readouterr().out)["kind"] == "ImportLabelsError"
    assert not pairs.exists()


# ---- the folder: merged rows, the command, the card ---------------------------------------------------------

def test_imports_merge_into_the_pairs_folder_and_the_card_says_where_they_came_from(images, tmp_path, capsys):
    from germquant.import_labels import import_labels

    lab2, lab3 = _labels_2d(), _labels_3d()
    img2 = images.add("cells_001.nd2", _pixels(lab2, 1)[:, None], SP2, two_d=True)
    img3 = images.add("stack_001.nd2", _pixels(lab3, 2), SP3, two_d=False)
    pairs = tmp_path / "pairs"
    m2 = _write_labels(tmp_path / "l" / "m2.tif", "tiff", lab2)
    m3 = _write_labels(tmp_path / "l" / "m3.npz", "npz", lab3)
    argv = ["import-labels", str(m2), "--format", "tiff", "--image", str(img2), "--pairs", str(pairs), "--json"]
    assert cli.main(argv + ["--channel", "1"]) == 0
    rec = json.loads(capsys.readouterr().out)
    assert rec["channel_index"] == 1 and rec["channel_name"] == "foci" and rec["n_pairs"] == 1
    assert np.array_equal(tifffile.imread(pairs / rec["pairs"][0]["image"]),
                          _u16(images.data("cells_001.nd2")[1, 0]))
    # the same image again replaces its rows; another image's rows join them
    assert cli.main(argv) == 0
    assert json.loads(capsys.readouterr().out)["replaced"] == 1
    res3 = import_labels(m3, "cellacdc-npz", img3, pairs, tool="micro-sam")
    man = _manifest(pairs)
    assert list(man["image_id"]).count("cells_001") == 1 and res3["replaced"] == 0
    assert set(man["source"]) == {"import:tiff", "import:micro-sam"}
    assert len(man) == 1 + res3["n_pairs"]
    # the text summary of the command
    assert cli.main(argv[:-1]) == 0
    out = capsys.readouterr().out
    assert "cells_001" in out and "import:tiff" in out and "manifest.csv" in out
    # the card: 2D and 3D pairs, the pixel sizes of both images, no pack (the labels came from no run)
    pairs_list = TR.load_pairs(pairs)
    card = TR._card_fields(pairs, pairs_list)
    assert card["dims"] == "2d+3d" and card["pixel_um"] is None and card["pixel_um_values"] == [0.2, 0.65]
    assert card["pack"] is None and "imported" in card["pack_note"] and "import:micro-sam" in card["pack_note"]
    assert card["segmented_role"] is None and card["object_kind"] == "nucleus"
    only2d = [p for p in pairs_list if p["image_id"] == "cells_001"]
    assert TR._card_fields(pairs, only2d)["dims"] == "2d"
    # an export into the same folder replaces the rows of the image it exports (stack_001, run in 3D) and keeps
    # the imported rows of the other image with their source and tile columns
    from germquant.config import load_config

    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("spots.enabled", False)
    pipeline.process_image(img3, cfg, tmp_path / "results" / "stack_001", prov=None)
    exported = TR.export_training_pairs(tmp_path / "results", pairs)
    assert list(exported.columns) == TR.MANIFEST_COLS + ["source", "tile_y0", "tile_x0"]
    man = _manifest(pairs)
    kept = man[man["image_id"] == "cells_001"]
    assert len(kept) == 1 and kept["source"].iloc[0] == "import:tiff" and int(kept["tile_y0"].iloc[0]) == 0
    new = man[man["image_id"] == "stack_001"]
    assert len(new) and new["source"].isna().all() and new["tile_y0"].isna().all()
    assert set(new["orientation"]) <= {"xy", "xz", "yz"} and set(new["channel_role"]) == {"dna"}
    raw = pd.read_csv(pairs / "manifest.csv", dtype=str, keep_default_na=False)
    assert raw.loc[raw["image_id"] == "cells_001", ["tile_y0", "tile_x0"]].values.tolist() == [["0", "0"]]
    assert set(raw.loc[raw["image_id"] == "stack_001", "tile_y0"]) == {""}             # whole offsets stay whole
    # the help says what the command does and points at the training guide
    capsys.readouterr()
    with pytest.raises(SystemExit):
        cli.main(["import-labels", "--help"])
    shown = capsys.readouterr().out
    assert "import:<tool>" in shown and "docs/TRAINING.md" in shown and "--allow-pickle" in shown


# ---- D28: the pipeline (torch) before the image is read ----------------------------------------------------

def test_import_labels_imports_the_pipeline_before_it_reads_the_image(tmp_path):
    """In a fresh interpreter, import_labels imports germquant.pipeline, and with it torch, before it reads an
    OME-TIFF through bioio, whose dask import loads CuPy (when installed) before torch otherwise (R5).
    Importing germquant.import_labels itself loads none of them."""
    pytest.importorskip("torch")
    pytest.importorskip("bioio_ome_tiff")
    import planar_helpers as PH

    img = PH.write_ome_cyx(tmp_path / "plane_001.ome.tif")
    info = PH.planted()["nuclei"]
    lab = _disks((PH.SIZE, PH.SIZE), [(n["y"], n["x"]) for n in info], 6, range(1, len(info) + 1))
    src = _write_labels(tmp_path / "plane_001_masks.tif", "tiff", lab)
    code = ("import json, sys\n"
            "import germquant.import_labels as IL\n"
            "light = [m for m in ('torch', 'germquant.pipeline', 'bioio', 'dask', 'cupy') if m in sys.modules]\n"
            f"res = IL.import_labels({str(src)!r}, 'tiff', {str(img)!r}, {str(tmp_path / 'pairs')!r})\n"
            "order = list(sys.modules)\n"
            "pos = {m: order.index(m) for m in ('torch', 'germquant.pipeline', 'bioio', 'dask', 'dask.array', 'cupy')\n"
            "       if m in sys.modules}\n"
            "print(json.dumps({'light': light, 'pos': pos, 'n': res['n_pairs'], 'dims': res['dims']}))\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300,
                       env=subprocess_env(), cwd=tmp_path)
    assert r.returncode == 0, r.stderr[-3000:]
    out = json.loads(r.stdout.strip().splitlines()[-1])
    assert out["light"] == [] and out["n"] == 1 and out["dims"] == "2d"
    pos = out["pos"]
    assert "torch" in pos and "germquant.pipeline" in pos and "bioio" in pos, pos
    readers = [m for m in ("bioio", "dask", "dask.array", "cupy") if m in pos]
    assert all(pos["torch"] < pos[m] and pos["germquant.pipeline"] < pos[m] for m in readers), pos
