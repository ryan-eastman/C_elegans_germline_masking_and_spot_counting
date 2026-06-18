"""QC flagging — unattended runs must surface problems, not hide them."""
from __future__ import annotations


def qc_flags(*, n_nuclei: int, channel_flags: list[str], axis_flags: list[str],
             spots_found: bool, spots_enabled: bool = True) -> tuple[bool, list[str]]:
    flags: list[str] = []
    flags += channel_flags + axis_flags

    if n_nuclei == 0:
        flags.append("qc:NO_NUCLEI")
    elif n_nuclei < 10:
        flags.append(f"qc:few_nuclei_{n_nuclei}")

    if not spots_enabled:
        flags.append("spots:disabled_segmentation_only")  # ran by choice (--no-spots), not a failure
    elif not spots_found:
        flags.append("qc:no_spots_detected")

    hard_fail = any("MISSING_REQUIRED" in f or "NO_NUCLEI" in f for f in flags)
    return (not hard_fail), flags
