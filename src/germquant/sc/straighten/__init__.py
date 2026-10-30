"""SC ribbon straightening and quantification (``germquant sc-straighten``; docs/SC_STRAIGHTEN.md).

A post-run command on finished fixed-image results: every synaptonemal-complex ribbon (the central_element
channel, SYP-3 or GFP::SYP-2) is traced in 3D inside its nucleus territory, its centreline refined and
fitted, its ends placed by an anchored half-maximum rule, and it is resampled along the centreline into a
straight strip (Fiji's Straighten, in 3D). Measurements along it: length, intensity profile, width, gaps
and pieces, the positions of foci (RAD-51 or COSA-1 spots) along the straightened ribbon, QC flags and an
abstaining strand class; a per-nucleus tier (totals, continuity, on-strand SYP fraction) that does not need a
perfect decomposition. Everything it writes is new files; nothing existing is touched.

Experimental: every image carries the flag ``sc_straighten:uncalibrated`` until the pre-registered protocol
(docs/SC_VALIDATION_PROTOCOL.md) has scored the method.
"""
from __future__ import annotations

VERSION = "2.0"

__all__ = ["VERSION", "run_results", "straighten_image"]


def __getattr__(name):
    if name == "straighten_image":
        from . import pipeline as _p

        return _p.straighten_image
    if name == "run_results":
        from . import runner as _r

        return _r.run_results
    raise AttributeError(name)
