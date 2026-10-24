"""sc-straighten inputs: the channel-subset dask construction (the same as germquant.io.read_stack), labels read
page by page into uint16, discovery by image_id with ids that end in an underscore, and on the workstation the
equality with read_stack on a development image (skipped when the file is absent)."""
from __future__ import annotations

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


@pytest.mark.skipif(not HERM000.is_file(), reason="development image on the lab workstation only")
def test_channel_subset_equals_read_stack_on_a_development_image():
    from germquant.io import read_stack

    st = read_stack(HERM000, z_range=(30, 33))
    got = inputs.nd2_role_arrays(HERM000, {"central_element": 1, "foci": 2}, z_range=[30, 33])
    assert np.array_equal(got["central_element"], st.data[1]) and np.array_equal(got["foci"], st.data[2])
    assert got["central_element"].dtype == st.data.dtype
    planes = inputs.nd2_role_planes(HERM000, {"foci": 2}, z_range=[30, 33])
    assert planes.shape[1] == 3 and np.array_equal(planes("foci", 1), st.data[2][1])
