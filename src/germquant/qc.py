"""QC flagging — unattended runs must surface problems, not hide them."""
from __future__ import annotations


def qc_flags(*, n_nuclei: int, channel_flags: list[str], axis_flags: list[str],
             zone_flags: list[str], sc_traced: bool, mean_sc_len, foci_found: bool) -> tuple[bool, list[str]]:
    flags: list[str] = []
    flags += channel_flags + axis_flags + zone_flags

    if n_nuclei == 0:
        flags.append("qc:NO_NUCLEI")
    elif n_nuclei < 10:
        flags.append(f"qc:few_nuclei_{n_nuclei}")

    if not sc_traced:
        flags.append("qc:no_sc_tracing")
    elif mean_sc_len is not None and mean_sc_len == mean_sc_len:  # not NaN
        # implausible total SC length per pachytene nucleus (WT ~ tens of µm)
        if mean_sc_len > 200:
            flags.append(f"qc:implausible_sc_length_{mean_sc_len:.0f}um")

    if not foci_found:
        flags.append("qc:no_foci_detected")

    hard_fail = any(
        "MISSING_REQUIRED" in f or "NO_NUCLEI" in f for f in flags
    )
    return (not hard_fail), flags
