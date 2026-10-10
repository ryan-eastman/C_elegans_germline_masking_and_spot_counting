"""Time-lapse quantification (docs/DESIGN_timelapse.md, docs/TIMELAPSE.md).

The per-timepoint loop lives in ``series`` (process_series: every enabled stage on each timepoint of one
time-lapse file, the series record, the series_summary table and the reference volumes for registration);
the drift registration in ``register`` (register_series, ``germquant register``: the register_shifts table and
the registration report). Tracking and the biologist validation come in later packages and write the other
schema.SERIES_TABLES.
"""
from .register import ReferenceVolumesMissing, register_series
from .series import SERIES_RECORD_SUFFIX, parse_timepoints, process_series, read_series_record

__all__ = ["SERIES_RECORD_SUFFIX", "ReferenceVolumesMissing", "parse_timepoints", "process_series",
           "read_series_record", "register_series"]
