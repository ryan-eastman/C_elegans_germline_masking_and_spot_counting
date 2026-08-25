"""Canonical self-verifying channel loader for ccw77 (PERMANENT copy in coloc_analysis so it survives
temp-scratchpad cleanup). Channel identity verified 3 ways, fixed BY INDEX:
    0 = 405/438 DAPI | 1 = 477/511 PGL-1(GFP) | 2 = 545/595 SYP-3(mCherry) | 3 = 640/681 lamin
SP = (dz,dy,dx) = (0.2, 0.1083, 0.1083) um. load() -> role -> (Z,Y,X) array."""
import sys

import numpy as np

sys.path.insert(0, r"C:\Users\ryane\C_elegans_germline_masking_and_spot_counting\src")
from germquant.io import read_stack

SP = np.array([0.2, 0.1083, 0.1083])
IDX = {"dapi": 0, "pgl": 1, "syp": 2, "lamin": 3}


# acceptable nd2 channel names per index (laser lines); anything else = acquisition order changed
EXPECT = {0: ("405", "DAPI"), 1: ("477", "488"), 2: ("545", "555", "561"), 3: ("640", "647")}


def verify_names(names, path):
    """raise if the nd2 channel order does not match the assumed DAPI/PGL/SYP/lamin index layout.
    This guard existed in the original loader, was lost in a rebuild, and is now restored: a silent
    channel swap here would put granule segmentation on the SYP channel and invert the biology."""
    if len(names) < 4:
        raise ValueError(f"chload: expected >=4 channels, got {names} in {path}")
    for i, ok in EXPECT.items():
        nm = str(names[i])
        if not any(k.lower() in nm.lower() for k in ok):
            raise ValueError(f"chload: channel {i} is {nm!r}, expected one of {ok} in {path}. "
                             "Acquisition channel order changed: do NOT analyze until IDX is re-verified.")


def load(path, verify_spacing=True):
    p = str(path).replace(chr(92), "/")
    st = read_stack(p)
    data = st.data  # (C, Z, Y, X)
    verify_names(st.channel_names, p)
    if verify_spacing:
        if not st.spacing_ok:
            raise ValueError(f"chload: voxel size unreadable in {p}; cannot trust micron-scale analysis")
        dz, dy, dx = st.spacing
        if not (abs(dz - SP[0]) < 0.03 and abs(dy - SP[1]) < 0.01 and abs(dx - SP[2]) < 0.01):
            raise ValueError(f"chload: spacing {st.spacing} != expected {tuple(SP)} for {p}")
    return {role: data[i] for role, i in IDX.items()}


# ---------------------------------------------------------------------------------------------
# image-id parsing + audit exclusions (single source of truth: coloc_analysis/exclusions.json)
import json as _json
import os as _os
import re as _re

_EXCL_PATH = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "exclusions.json")


def parse_iid(iid):
    """-> dict(batch, sex, treat, short). short = '<treat>_<sex>_<num>', e.g. 'HS_male_07'."""
    m = _re.search(r"(noHS|HS)_(male|herm)_?(\d+)", iid)
    short = f"{m.group(1)}_{m.group(2)}_{m.group(3)}" if m else iid
    batch = iid[:8] if _re.match(r"\d{8}_", iid) else "?"
    return {"batch": batch, "sex": m.group(2) if m else "?", "treat": m.group(1) if m else "?", "short": short}


def is_excluded(iid):
    """True if this gonad is in the audit exclusion list (exclusions.json)."""
    ex = _json.load(open(_EXCL_PATH, encoding="utf-8"))
    p = parse_iid(iid)
    return p["batch"] == ex["excluded_batch"] and p["short"] in set(ex["excluded_short_ids"])
