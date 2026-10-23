"""Per-nucleus crops, the nucleus territory and the label checks.

One ``ndi.find_objects`` pass; each box is padded by ``crop_pad_um`` (z, y, x) and the crops are copies (a few
hundred kB each), so a worker process never needs the whole image. The territory is the set of voxels within
``territory_reach_um`` of the nucleus's own label (Euclidean distance in microns) that lie in no other label:
the median 3 percent of ridge signal that the DAPI label misses at the nuclear periphery stays reachable.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage as ndi


@dataclass
class Crop:
    nucleus_id: int
    org: tuple                    # whole-image voxel index of crop[0, 0, 0]
    ce: np.ndarray                # central_element crop (raw counts, the image dtype)
    labels: np.ndarray            # label crop (uint16 or int32)
    roles: dict = field(default_factory=dict)   # other roles' crops (profiles, strips)
    touches_z: bool = False
    touches_xy: bool = False
    image_shape: tuple = ()

    @property
    def shape(self) -> tuple:
        return tuple(self.ce.shape)


def boxes(labels: np.ndarray) -> dict[int, tuple]:
    """{label: slices} for every label present (``ndi.find_objects``)."""
    out = {}
    for i, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is not None:
            out[i] = sl
    return out


def pad_vox(spacing, pad_um) -> tuple[int, int, int]:
    return tuple(int(np.ceil(float(p) / float(s) - 1e-9)) for p, s in zip(pad_um, spacing))


def make_crop(nid: int, box: tuple, ce: np.ndarray, labels: np.ndarray, spacing, pad_um,
              roles: dict | None = None) -> Crop:
    shape = labels.shape
    pv = pad_vox(spacing, pad_um)
    sl = tuple(slice(max(0, b.start - p), min(n, b.stop + p)) for b, p, n in zip(box, pv, shape))
    org = tuple(int(s.start) for s in sl)
    touches_z = box[0].start == 0 or box[0].stop == shape[0]
    touches_xy = box[1].start == 0 or box[2].start == 0 or box[1].stop == shape[1] or box[2].stop == shape[2]
    rc = {k: np.array(v[sl]) for k, v in (roles or {}).items() if v is not None}
    return Crop(nucleus_id=int(nid), org=org, ce=np.array(ce[sl]), labels=np.array(labels[sl]), roles=rc,
                touches_z=bool(touches_z), touches_xy=bool(touches_xy), image_shape=tuple(shape))


def territory(lab: np.ndarray, nid: int, spacing, reach_um: float) -> np.ndarray:
    own = lab == nid
    if not own.any():
        return own
    dist = ndi.distance_transform_edt(~own, sampling=spacing)
    return (dist <= reach_um) & ((lab == 0) | own)


def other_label_distance(lab: np.ndarray, nid: int, spacing) -> np.ndarray:
    """Distance (um) to the nearest voxel of another nonzero label (inf when there is none in the crop)."""
    other = (lab != 0) & (lab != nid)
    if not other.any():
        return np.full(lab.shape, np.inf, np.float32)
    return ndi.distance_transform_edt(~other, sampling=spacing).astype(np.float32)
