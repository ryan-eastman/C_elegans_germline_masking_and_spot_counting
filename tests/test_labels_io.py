"""germquant.io.labels: the readers of segmentation.nuclei.labels_from (docs/DESIGN_system_packs.md 4.11, decisions
D13 and D14 of the Train B plan; docs/LABELS.md). No pipeline run here; tests/test_labels_run.py runs them.

* TIFF, OME-TIFF and Cell-ACDC segm.npz files are read as 2D, 3D and with a T axis (labels_frame picks the frame),
  always to int32 (Z, Y, X) relabelled 1..N in the order of the file's ids;
* xy_stride and z_range cut the labels as the image readers cut the stack; a 2D file is copied through z;
* npz files are read with allow_pickle=False, so a pickled npz is refused; a Cellpose _seg.npy path is refused on
  the run path with the import-labels hint, and read_seg_npy reads one only with allow_pickle;
* labels_path formats the template ({image_dir}, {image_stem}, {image_id}); recipe_keys checks the two keys, the
  pack lint allow-list carries labels_frame, and check_planar skips the segmentation backend in labels runs.
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pytest
from labels_helpers import (
    write_npz,
    write_ome,
    write_pickled_npz,
    write_seg_npy,
    write_tiff,
)

from germquant import planar, recipe_keys
from germquant.config import load_config
from germquant.io import labels as L

ROOT = Path(__file__).resolve().parents[1]
NZ, NY, NX = 6, 12, 16


def _lab3d() -> np.ndarray:
    """(NZ, NY, NX) labels with ids that are not consecutive (5, 9, 40), one object cut off by z_range (2, 5)."""
    lab = np.zeros((NZ, NY, NX), np.int32)
    lab[0:2, 1:4, 1:5] = 40          # only in planes 0 and 1
    lab[1:5, 5:9, 2:7] = 9
    lab[2:6, 2:6, 9:14] = 5
    return lab


def _relabelled(lab: np.ndarray) -> np.ndarray:
    """The expected relabel: the ids present, in increasing order, become 1..N."""
    out = np.zeros(lab.shape, np.int32)
    for new, old in enumerate(sorted(set(np.unique(lab)) - {0}), 1):
        out[lab == old] = new
    return out


@pytest.mark.parametrize("kind", ["tiff", "ome", "npz2d", "npz3d", "npzT"])
def test_reads_labels(tmp_path, kind):
    lab3 = _lab3d()
    lab2 = lab3[2]
    series = np.stack([lab3, np.where(lab3 > 0, lab3 + 100, 0)])          # (T, Z, Y, X): frame 1 differs
    if kind == "tiff":
        cases = [(write_tiff(tmp_path / "a2.tif", lab2), 0, lab2, NZ),
                 (write_tiff(tmp_path / "a3.tif", lab3), 0, lab3, NZ),
                 (write_tiff(tmp_path / "a3u.tiff", lab3.astype(np.uint16)), 0, lab3, NZ)]
        # tifffile stores a stack of 3 or 4 planes as one page of separate samples (axes SYX): still planes
        small = lab3[:3]
        cases.append((write_tiff(tmp_path / "a3s.tif", small), 0, small, 3))
    elif kind == "ome":
        cases = [(write_ome(tmp_path / "b2.ome.tif", lab2, "YX"), 0, lab2, NZ),
                 (write_ome(tmp_path / "b3.ome.tif", lab3, "ZYX"), 0, lab3, NZ),
                 (write_ome(tmp_path / "bt.ome.tiff", series, "TZYX"), 1, series[1], NZ),
                 (write_ome(tmp_path / "bty.ome.tif", series[:, 2], "TYX"), 1, series[1, 2], NZ)]
    elif kind == "npz2d":
        cases = [(write_npz(tmp_path / "c2_segm.npz", lab2), 0, lab2, NZ)]
    elif kind == "npz3d":
        cases = [(write_npz(tmp_path / "c3_segm.npz", lab3), 0, lab3, NZ)]
    else:
        # Cell-ACDC time-lapse layouts: (T, Z, Y, X), and (T, Y, X) for a series of 2D labels
        cases = [(write_npz(tmp_path / "ct_segm.npz", series), 1, series[1], NZ),
                 (write_npz(tmp_path / "ct2_segm.npz", series[:, 2]), 0, series[0, 2], NZ),
                 (write_npz(tmp_path / "ct3_segm.npz", series[:, 2]), 1, series[1, 2], NZ)]
    for path, frame, want, nz in cases:
        got = L.read_labels(path, nz=nz, frame=frame)
        assert got.dtype == np.int32 and got.shape == (nz, NY, NX), (path.name, got.shape)
        expect = _relabelled(want if want.ndim == 3 else np.repeat(want[None], nz, axis=0))
        np.testing.assert_array_equal(got, expect, err_msg=path.name)
        assert int(got.max()) == len(set(np.unique(want)) - {0})


def test_stride_and_z_range(tmp_path):
    lab3 = _lab3d()
    p3 = write_tiff(tmp_path / "s3.tif", lab3)
    got = L.read_labels(p3, nz=NZ, xy_stride=2, z_range=(2, 5))
    np.testing.assert_array_equal(got, _relabelled(lab3[2:5, ::2, ::2]))
    assert int(got.max()) == 2                      # id 40 lives in planes 0 and 1 only: cut away, then 1..N
    # a 2D file: the stride applies and the plane is copied through the planes z_range keeps
    p2 = write_npz(tmp_path / "s2_segm.npz", lab3[3])
    got2 = L.read_labels(p2, nz=NZ, xy_stride=3, z_range=(1, 5))
    assert got2.shape == (4, 4, 6)
    np.testing.assert_array_equal(got2, np.repeat(_relabelled(lab3[3, ::3, ::3])[None], 4, axis=0))
    # z_range past the last plane keeps what exists, as the image readers' slice does
    assert L.read_labels(p3, nz=NZ, z_range=(4, 99)).shape == (2, NY, NX)
    assert L.read_labels(p2, nz=NZ, z_range=(4, 99)).shape == (2, NY, NX)


def test_2d_labels_copied_through_z_and_relabelled(tmp_path):
    lab2 = np.zeros((NY, NX), np.uint16)
    lab2[1:4, 1:4], lab2[6:9, 2:6], lab2[2:5, 10:14] = 42, 7, 3000
    for path in (write_tiff(tmp_path / "p.tif", lab2), write_npz(tmp_path / "p_segm.npz", lab2),
                 write_ome(tmp_path / "p.ome.tif", lab2[None], "ZYX")):   # one plane: a 2D file too
        got = L.read_labels(path, nz=9)
        assert got.shape == (9, NY, NX) and got.dtype == np.int32
        assert all(np.array_equal(got[z], got[0]) for z in range(9))
        assert set(np.unique(got)) == {0, 1, 2, 3}
        assert got[0, 7, 3] == 1 and got[0, 2, 2] == 2 and got[0, 3, 12] == 3      # 7, 42, 3000 in order
    # no size gate: a one-pixel object stays
    tiny = np.zeros((NY, NX), np.int32)
    tiny[5, 5] = 2
    assert int(L.read_labels(write_tiff(tmp_path / "tiny.tif", tiny), nz=1).max()) == 1
    # an all-zero file is a labels file with no object
    assert not L.read_labels(write_tiff(tmp_path / "zero.tif", np.zeros((NY, NX), np.uint8)), nz=2).any()
    # integral floats and booleans are labels; fractions and negative ids are not
    assert L.read_labels(write_npz(tmp_path / "f_segm.npz", lab2.astype(np.float32)), nz=1).max() == 3
    assert L.read_labels(write_npz(tmp_path / "b_segm.npz", lab2 > 0), nz=1).max() == 1
    for name, bad in (("frac", lab2.astype(np.float32) + 0.5), ("neg", lab2.astype(np.int32) * -1)):
        with pytest.raises(L.LabelsError, match="label ids"):
            L.read_labels(write_npz(tmp_path / f"{name}_segm.npz", bad), nz=1)


def test_the_file_must_match_the_images_planes(tmp_path):
    lab3 = _lab3d()
    p = write_ome(tmp_path / "z.ome.tif", lab3, "ZYX")
    with pytest.raises(L.LabelsError, match=r"6 planes.*the image has 8"):
        L.read_labels(p, nz=8)
    # a segm.npz or plain TIFF whose first axis is not the image's planes is a series of 2D labels (T, Y, X)
    for path in (write_npz(tmp_path / "amb_segm.npz", lab3), write_tiff(tmp_path / "amb.tif", lab3)):
        np.testing.assert_array_equal(L.read_labels(path, nz=6), _relabelled(lab3))           # 6 = the image's Z
        np.testing.assert_array_equal(L.read_labels(path, nz=8, frame=2),
                                      np.repeat(_relabelled(lab3[2])[None], 8, axis=0))      # frame 2 of 6
        np.testing.assert_array_equal(L.read_labels(path, nz=1, frame=4), _relabelled(lab3[4])[None])
    # what a label image cannot be
    with pytest.raises(L.LabelsError, match="channel"):
        L.read_labels(write_ome(tmp_path / "c.ome.tif", np.stack([lab3, lab3]), "CZYX"), nz=6)
    with pytest.raises(L.LabelsError, match="colour"):
        import tifffile

        rgb = tmp_path / "rgb.tif"
        tifffile.imwrite(str(rgb), np.zeros((NY, NX, 3), np.uint8), photometric="rgb")
        L.read_labels(rgb, nz=1)
    multi = tmp_path / "multi.tif"
    import tifffile

    with tifffile.TiffWriter(str(multi)) as w:
        w.write(lab3[0])
        w.write(np.zeros((5, 5), np.int32))
    with pytest.raises(L.LabelsError, match="2 images"):
        L.read_labels(multi, nz=1)
    with pytest.raises(L.LabelsError, match="axes"):
        L.read_labels(write_npz(tmp_path / "five_segm.npz", np.zeros((2, 2, 2, NY, NX), np.int32)), nz=2)


def test_labels_frame_picks_the_frame_of_a_t_axis(tmp_path):
    """D13: segmentation.nuclei.labels_frame names the frame of a T axis (0 by default)."""
    lab3 = _lab3d()
    frames = np.stack([lab3, np.roll(lab3, 3, axis=2), np.roll(lab3, 6, axis=2)])
    for path in (write_npz(tmp_path / "t_segm.npz", frames), write_ome(tmp_path / "t.ome.tif", frames, "TZYX")):
        for f in range(3):
            np.testing.assert_array_equal(L.read_labels(path, nz=NZ, frame=f), _relabelled(frames[f]))
        with pytest.raises(L.LabelsError, match=r"labels_frame 3 .* 3 frames \(0 to 2\)"):
            L.read_labels(path, nz=NZ, frame=3)
    # a file without a T axis has frame 0 only
    no_t = write_ome(tmp_path / "zonly.ome.tif", lab3, "ZYX")
    with pytest.raises(L.LabelsError, match="labels_frame 1 .* no T axis"):
        L.read_labels(no_t, nz=NZ, frame=1)
    # a T axis of one frame (tifffile drops it from an OME-TIFF's axes): frame 0 reads it, frame 1 does not exist
    for one in (write_ome(tmp_path / "t1.ome.tif", lab3[None], "TZYX"), write_npz(tmp_path / "t1_segm.npz", lab3[None])):
        np.testing.assert_array_equal(L.read_labels(one, nz=NZ), _relabelled(lab3))
        with pytest.raises(L.LabelsError, match=r"labels_frame 1 .*(has 1 frame \(0 to 0\)|has no T axis)"):
            L.read_labels(one, nz=NZ, frame=1)
    # the key: read only in labels runs, 0 by default
    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    assert L.labels_frame(cfg) == 0
    cfg.set("segmentation.nuclei.labels_frame", 2)
    assert L.labels_frame(cfg) == 2


def test_pickled_npz_refused(tmp_path):
    path = write_pickled_npz(tmp_path / "obj_segm.npz")
    with pytest.raises(L.LabelsError, match="allow_pickle=False"):
        L.read_labels(path, nz=1)
    # a pickle saved under an .npz name is refused too, never unpickled
    raw = tmp_path / "raw_segm.npz"
    np.save(str(tmp_path / "raw.npy"), np.array([{"a": 1}], dtype=object), allow_pickle=True)
    raw.write_bytes((tmp_path / "raw.npy").read_bytes())
    with pytest.raises(L.LabelsError, match="not an npz archive|allow_pickle=False"):
        L.read_labels(raw, nz=1)
    empty = tmp_path / "empty_segm.npz"
    np.savez(str(empty))
    with pytest.raises(L.LabelsError, match="no array"):
        L.read_labels(empty, nz=1)


def test_seg_npy_refused_with_import_labels_hint(tmp_path):
    masks = np.zeros((NY, NX), np.uint16)
    masks[2:6, 3:8] = 4
    seg = write_seg_npy(tmp_path / "img_seg.npy", masks)
    for path in (seg, tmp_path / "plain.npy"):
        with pytest.raises(L.LabelsError) as err:
            L.read_labels(path, nz=1)
        msg = str(err.value)
        assert "import-labels" in msg and "--allow-pickle" in msg and "never reads" in msg, msg
    with pytest.raises(L.LabelsError, match="import-labels"):
        L.labels_format(seg)
    for bad in ("x.png", "x.h5", "x.zarr", "x"):
        with pytest.raises(L.LabelsError, match=r"\.ome\.tif.*\.npz"):
            L.labels_format(tmp_path / bad)
    assert [L.labels_format(Path(n)) for n in ("a.TIF", "a.tiff", "a.ome.tif", "a.ome.tiff", "a_segm.npz")] == \
        ["tiff", "tiff", "tiff", "tiff", "npz"]
    # import-labels' reader (P5): only with allow_pickle, and the masks come back as they are
    with pytest.raises(L.LabelsError, match="allow_pickle"):
        L.read_seg_npy(seg, allow_pickle=False)
    np.testing.assert_array_equal(L.read_seg_npy(seg, allow_pickle=True), masks)
    plain = tmp_path / "plain_masks.npy"
    np.save(str(plain), masks)                       # not a Cellpose dict
    with pytest.raises(L.LabelsError, match="masks"):
        L.read_seg_npy(plain, allow_pickle=True)


def test_labels_path_template(tmp_path):
    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    image = tmp_path / "data" / "20260101_cells_001.ome.tif"
    assert L.labels_path(cfg, image) is None and L.labels_template(cfg) is None      # unset: not a labels run
    cfg.set("segmentation.nuclei.labels_from", "{image_dir}/masks/{image_stem}_masks.tif")
    assert L.labels_path(cfg, image) == tmp_path / "data" / "masks" / "20260101_cells_001_masks.tif"
    cfg.set("segmentation.nuclei.labels_from", "{image_dir}/{image_id}_segm.npz")
    assert L.labels_path(cfg, tmp_path / "x" / "g_HS_male_7.nd2") == tmp_path / "x" / "g_HS_male_7_segm.npz"
    # a relative image path: {image_dir} is its absolute folder
    rel = L.labels_path(cfg, Path("some") / "g.nd2")
    assert rel.is_absolute() and rel == (Path.cwd() / "some" / "g_segm.npz")
    # a template without {image_dir} that is relative resolves against the recipe's folder, as io.channel_map does
    cfg.set("segmentation.nuclei.labels_from", "truth/{image_id}.tif")
    assert L.labels_path(cfg, image) == cfg.base_dir / "truth" / "20260101_cells_001.tif"
    absolute = tmp_path / "abs" / "{image_stem}.tif"
    cfg.set("segmentation.nuclei.labels_from", absolute.as_posix())
    assert L.labels_path(cfg, image) == tmp_path / "abs" / "20260101_cells_001.tif"
    cfg.set("segmentation.nuclei.labels_from", "{image_dir}/{image_name}.tif")
    with pytest.raises(L.LabelsError, match="image_name"):
        L.labels_path(cfg, image)


def test_labels_sha256_is_the_file_digest(tmp_path):
    path = write_tiff(tmp_path / "d.tif", _lab3d())
    assert L.labels_sha256(path) == hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(FileNotFoundError):
        L.labels_sha256(tmp_path / "missing.tif")


BAD_KEYS = [
    ("segmentation.nuclei.labels_from", 3), ("segmentation.nuclei.labels_from", True),
    ("segmentation.nuclei.labels_from", ""), ("segmentation.nuclei.labels_from", "   "),
    ("segmentation.nuclei.labels_from", ["a.tif"]), ("segmentation.nuclei.labels_from", {"image_dir": 1}),
    ("segmentation.nuclei.labels_from", "{image_dir}/{stem}.tif"), ("segmentation.nuclei.labels_from", "{0}.tif"),
    ("segmentation.nuclei.labels_from", "{image_dir.parent}/x.tif"), ("segmentation.nuclei.labels_from", "{image_id"),
    ("segmentation.nuclei.labels_from", "{image_id:>9}.tif"), ("segmentation.nuclei.labels_from", "{image_id!r}.tif"),
    ("segmentation.nuclei.labels_frame", -1), ("segmentation.nuclei.labels_frame", 1.0),
    ("segmentation.nuclei.labels_frame", True), ("segmentation.nuclei.labels_frame", "2"),
    ("segmentation.nuclei.labels_frame", np.int64(1)),
]


@pytest.mark.parametrize("key,value", BAD_KEYS)
def test_recipe_keys_refuse_a_malformed_labels_key(key, value):
    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    cfg._data.setdefault("segmentation", {}).setdefault("nuclei", {})[key.rsplit(".", 1)[1]] = value
    (problem,) = recipe_keys.problems(cfg)
    assert problem.startswith(f"{key} must be"), problem
    with pytest.raises(ValueError, match=re.escape(key)):
        recipe_keys.check(cfg)


@pytest.mark.parametrize("key,value", [
    ("segmentation.nuclei.labels_from", "{image_dir}/{image_stem}_masks.tif"),
    ("segmentation.nuclei.labels_from", "D:/truth/{image_id}_segm.npz"),
    ("segmentation.nuclei.labels_from", "one_file_for_every_image.tif"),
    ("segmentation.nuclei.labels_from", "{image_dir}/{{literal}}_{image_id}.tif"),
    ("segmentation.nuclei.labels_frame", 0), ("segmentation.nuclei.labels_frame", 12),
    ("segmentation.nuclei.labels_from", None), ("segmentation.nuclei.labels_frame", None),
])
def test_recipe_keys_accept_a_valid_labels_key(key, value):
    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    cfg._data.setdefault("segmentation", {}).setdefault("nuclei", {})[key.rsplit(".", 1)[1]] = value
    assert recipe_keys.problems(cfg) == []
    assert key in recipe_keys.CHECKED_KEYS


def test_the_lint_allow_list_carries_both_keys():
    from germquant.packs import lint

    for key in ("segmentation.nuclei.labels_from", "segmentation.nuclei.labels_frame"):
        assert key in lint.T5_KEYS and lint.key_readable(key, set())


def test_check_planar_skips_the_segmentation_backend_in_labels_runs():
    """D14: nothing is segmented in a labels run, so a backend without planar_ok is not refused there."""
    roles = {"dna": 0, "foci": 1, "granule": None}
    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    cfg.set("segmentation.nuclei.method", "cellpose")
    cfg.set("segmentation.nuclei.backend", "nuclei3d")
    with pytest.raises(planar.PlanarConfigError, match="'nuclei3d' does not declare planar_ok"):
        planar.check_planar(cfg, roles)
    cfg.set("segmentation.nuclei.labels_from", "{image_dir}/{image_stem}_masks.tif")
    planar.check_planar(cfg, roles)
    # the stage refusals still hold in a labels run
    cfg.set("axis.enabled", True)
    with pytest.raises(planar.PlanarConfigError, match="the axis stage"):
        planar.check_planar(cfg, roles)
