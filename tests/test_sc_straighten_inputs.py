"""sc-straighten inputs: the channel-subset dask construction (the same as germquant.io.read_stack), labels read
page by page into uint16, discovery by image_id with ids that end in an underscore, and on the workstation the
equality with read_stack on a development image (skipped when the file is absent)."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from germquant.sc.straighten import inputs

HERM000 = Path("C:/Users/ryane/n2_rad51_local/20251105_N2_noHS/20251105_N2_nohs_HERM _.nd2")


def test_canonical_dask_dims():
    da = pytest.importorskip("dask.array")
    raw = np.arange(2 * 4 * 3 * 5 * 6).reshape(2, 4, 3, 5, 6)          # T, Z, C, Y, X
    d = inputs.canonical_dask(da.from_array(raw), ["T", "Z", "C", "Y", "X"], xy_stride=2, z_range=[1, 3])
    out = np.asarray(d)
    assert out.shape == (3, 2, 3, 3)
    assert np.array_equal(out, raw[0].transpose(1, 0, 2, 3)[:, 1:3, ::2, ::2])
    d2 = inputs.canonical_dask(da.from_array(raw[0, 0]), ["C", "Y", "X"])  # no Z axis: one plane
    assert np.asarray(d2).shape == (3, 1, 5, 6)


def test_labels_page_by_page_into_uint16(tmp_path):
    tifffile = pytest.importorskip("tifffile")
    lab = np.zeros((4, 8, 9), np.int32)
    lab[1:3, 2:5, 2:6] = 7
    lab[3, 6:, 6:] = 70000
    p = tmp_path / "l.tif"
    tifffile.imwrite(p, lab, compression="zlib", ome=True, metadata={"axes": "ZYX"})
    small = lab.copy()
    small[small == 70000] = 9
    p2 = tmp_path / "s.tif"
    tifffile.imwrite(p2, small, compression="zlib", ome=True, metadata={"axes": "ZYX"})
    a = inputs.read_labels(p2, max_label=9)
    assert a.dtype == np.uint16 and np.array_equal(a, small)
    b = inputs.read_labels(p, max_label=70000)
    assert b.dtype == np.int32 and np.array_equal(b, lab)
    c = inputs.read_labels(p, max_label=9)                         # an understated maximum: promoted to int32
    assert np.array_equal(c, lab)


def test_discovery_by_image_id(tmp_path):
    d = tmp_path / "r" / "n2rad51_noHS_herm_000"
    d.mkdir(parents=True)
    iid = "20251105_N2_nohs_HERM _"
    pd.DataFrame({"nucleus_id": [1], "image_id": [iid]}).to_csv(d / f"{iid}__nuclei.csv", index=False)
    pd.DataFrame({"x": [1]}).to_csv(d / f"{iid}__sc_strand_nuclei.csv", index=False)   # never mistaken for one
    found = inputs.discover(tmp_path / "r")
    assert [f.image_id for f in found] == [iid]
    assert inputs.select(found, ["n2rad51_noHS_herm_000"]) == found
    assert inputs.select(found, [iid]) == found and inputs.select(found, ["other"]) == []


def test_an_image_sealed_by_its_name_is_refused_before_any_file_is_read(tmp_path, monkeypatch):
    # protocol changelog entry 27: discovery reads no nuclei.csv of an image sealed by its file name or folder name
    # (a test or reserved image stem or RAD-51 tag, without regard to case), and the runner refuses it before
    # reading any of its files; the files here are dummies that cannot be parsed
    from germquant.benchmark import sc_split
    from germquant.sc.straighten import params, runner

    g = next(x for x in sc_split.load_split()["gonads"] if x.get("role") == "test" and x.get("rad51_tag"))
    root = tmp_path / "r"
    for folder, name in ((g["rad51_tag"], g["image_stem"]), ("copies", g["image_stem"].upper()),
                         (g["rad51_tag"].upper(), "renamed")):
        (root / folder).mkdir(parents=True, exist_ok=True)
        (root / folder / f"{name}__nuclei.csv").write_bytes(b"\x00 not a table")
    dev = root / "n2rad51_noHS_herm_000"
    dev.mkdir()
    iid = "20251105_N2_nohs_HERM _"
    pd.DataFrame({"nucleus_id": [1], "image_id": [iid]}).to_csv(dev / f"{iid}__nuclei.csv", index=False)
    reads = []
    orig = pd.read_csv
    monkeypatch.setattr(pd, "read_csv", lambda path, *a, **k: reads.append(str(path)) or orig(path, *a, **k))
    found = inputs.discover(root)
    assert sorted(f.image_id for f in found) == sorted([g["image_stem"], g["image_stem"].upper(), "renamed", iid])
    assert len(reads) == 1 and reads[0].endswith(f"n2rad51_noHS_herm_000{os.sep}{iid}__nuclei.csv")
    sealed = [f for f in found if f.image_id != iid]
    assert all(inputs.sealed_by_name(f) == "test" for f in sealed)
    reads.clear()
    oc = runner.run_results(root, [f.image_id for f in sealed], params=params.resolve(), echo=lambda *_: None)
    # discovery reads the development image's table again; nothing of the sealed ones is read
    assert [o.code for o in oc] == ["sealed"] * 3 and all("n2rad51_noHS_herm_000" in r for r in reads)


@pytest.mark.skipif(not HERM000.is_file(), reason="development image on the lab workstation only")
def test_channel_subset_equals_read_stack_on_a_development_image():
    from germquant.io import read_stack

    st = read_stack(HERM000, z_range=(30, 33))
    got = inputs.nd2_role_arrays(HERM000, {"central_element": 1, "foci": 2}, z_range=[30, 33])
    assert np.array_equal(got["central_element"], st.data[1]) and np.array_equal(got["foci"], st.data[2])
    assert got["central_element"].dtype == st.data.dtype
    planes = inputs.nd2_role_planes(HERM000, {"foci": 2}, z_range=[30, 33])
    assert planes.shape[1] == 3 and np.array_equal(planes("foci", 1), st.data[2][1])
