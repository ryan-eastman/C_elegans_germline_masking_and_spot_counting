"""Light image header reads for projects (docs/PROJECTS.md, "The light entry"; PLAN.md decision D28).

Part 1a of unit P5-PRJ: the frozen signature, as a stub.
"""
from __future__ import annotations


def read_header(path) -> dict:
    """``{"sizes", "n_images", "channel_names", "spacing", "spacing_ok", "spacing_yx_um", "single_plane",
    "reader"}`` of one image, read without bioio for .nd2 and TIFF files."""
    raise NotImplementedError("germquant.project.headers.read_header is not built yet: unit P5-PRJ of tranche 5 "
                              "Train B builds it (C:/Users/ryane/germquant_golden/t5_train_b/PLAN.md, section 9.1)")
