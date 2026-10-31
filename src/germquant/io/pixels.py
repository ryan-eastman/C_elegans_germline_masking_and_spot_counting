"""The digest of an image's pixels, whatever file holds them (docs/TRAINING.md, step 3; tranche 5 package P5).

``pixels_sha256(array)`` is the sha256 of the ASCII text ``"<dtype>|<shape>|"`` (numpy's dtype.str, such as
``<u2``, and the sizes joined by commas, such as ``2,5,256,256``) followed by the array's bytes in C order. The
array is put in little-endian byte order first, so the same pixels give the same digest from a big-endian TIFF, on
any machine and from any container: an .nd2, a TIFF, an OME-TIFF or a copy in memory. A different dtype or shape
gives another digest.

germquant import-labels takes it of the image as germquant reads it (the (C, Z, Y, X) array at full resolution,
every plane, and for a time series the timepoint imported) and refuses an image whose digest, or the sha256 of
whose file, a committed benchmark split lists among its test images (benchmarks/ext/<dataset>/converted.json,
docs/DESIGN_system_packs.md section 6.6). A benchmark that writes converted.json takes the digest the same way.
"""
from __future__ import annotations

import hashlib

import numpy as np


def pixels_sha256(array) -> str:
    """sha256 of ``f"{dtype.str}|{shape joined by commas}|"`` (ASCII) and the C-order bytes of ``array``, in
    little-endian byte order (module docstring). ValueError for an array of Python objects, which has no bytes of
    its own to hash."""
    a = np.asarray(array)
    if a.dtype.hasobject:
        raise ValueError(f"pixels_sha256 takes an array of numbers, not of Python objects ({a.dtype})")
    if a.dtype.str.startswith(">"):
        a = a.astype(a.dtype.newbyteorder("<"))
    a = np.ascontiguousarray(a)
    h = hashlib.sha256(f"{a.dtype.str}|{','.join(str(int(n)) for n in a.shape)}|".encode("ascii"))
    h.update(memoryview(a.reshape(-1)).cast("B"))
    return h.hexdigest()
