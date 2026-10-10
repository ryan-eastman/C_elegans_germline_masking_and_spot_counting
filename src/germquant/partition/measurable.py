"""Partition measurability: did the SYP channel carry enough signal for a partition coefficient to mean
anything? (docs/PARTITION_QC.md)

The partition coefficient divides the SYP level in granule voxels by the SYP level in cytoplasm at the
same distance from the nuclear envelope. When the SYP channel holds no live fluorophore (a dead
mCherry::SYP-3 imaged at full laser power and a long exposure), or was recorded through the wrong
emission filter so that a far-red lamin stain leaks into it, both sides of that ratio are camera noise
plus whatever structured bleed sits on the envelope, and the coefficient settles on a floor above 1 that
neither the 180 degree rotation null nor the z-shift floor removes: the controls sample the same noise.
The numbers are then reproducible and meaningless. This module turns the SYP signal-to-noise that the
signal_qc stage writes to image_summary (``sig_central_element_snr``) into two image_summary columns
and, when a profile asks for it, one QC flag:

* ``partition_syp_snr``: that signal-to-noise; NaN when the field is absent, None or NaN (the signal_qc
  stage off, or the SYP role not resolved);
* ``partition_measurable``: ``partition_syp_snr >= partition.min_syp_snr`` when the gate is on
  (``partition.min_syp_snr`` above 0; the code default 0.0 keeps it off) and the ratio is known, else
  True;
* the flag ``partition:low_syp_snr`` when the gate is on and the image is not measurable. Every
  partition number is still written: the flag, not a blank, is the signal, so a reader can see what a
  dead channel produces and filter on the column.

On the default path (no ``partition.min_syp_snr`` in the config) only the two columns are added and no
flag fires, so the qc_flags text of every existing run is unchanged (docs/ROADMAP_modular_pipeline.md).
"""
from __future__ import annotations

import math
from collections.abc import Mapping

SYP_SNR_FIELD = "sig_central_element_snr"   # the signal_qc stage's image_summary field for the SYP role
LOW_SNR_FLAG = "partition:low_syp_snr"
MIN_SYP_SNR_DEFAULT = 0.0                   # code default of partition.min_syp_snr: the gate is off
MEASURABILITY_FIELDS = ("partition_syp_snr", "partition_measurable")


def syp_snr_from_summary(summary_fields: Mapping) -> float:
    """``sig_central_element_snr`` as a float; NaN when the field is absent, None, NaN or not a number."""
    try:
        value = float(summary_fields.get(SYP_SNR_FIELD))
    except (TypeError, ValueError, AttributeError):
        return float("nan")
    return value


def partition_measurability(summary_fields: Mapping, *,
                            min_syp_snr: float | None = MIN_SYP_SNR_DEFAULT) -> tuple[dict, list[str]]:
    """The two image_summary fields and the flags (empty, or ``[LOW_SNR_FLAG]``) for one image.

    ``summary_fields`` is the per-image summary dict process_image has accumulated by the time the
    partition stage has run. ``min_syp_snr`` is ``partition.min_syp_snr`` as the config holds it: None
    (a null in a profile) means off, like the code default; a value that is not a number raises, so a
    misspelt threshold is loud rather than silently off."""
    gate = MIN_SYP_SNR_DEFAULT if min_syp_snr is None else float(min_syp_snr)
    snr = syp_snr_from_summary(summary_fields)
    gate_on = gate > 0 and not math.isnan(snr)
    measurable = bool(snr >= gate) if gate_on else True
    fields = {"partition_syp_snr": snr, "partition_measurable": measurable}
    return fields, ([] if measurable else [LOW_SNR_FLAG])
