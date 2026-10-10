"""Time-lapse quantification (docs/DESIGN_timelapse.md, docs/TIMELAPSE.md).

The per-timepoint loop lives in ``series`` (process_series: every enabled stage on each timepoint of one
time-lapse file, the series record and the series_summary table). Registration, tracking and the biologist
validation come in later packages and write the other schema.SERIES_TABLES.
"""
from .series import SERIES_RECORD_SUFFIX, parse_timepoints, process_series, read_series_record

__all__ = ["SERIES_RECORD_SUFFIX", "parse_timepoints", "process_series", "read_series_record"]
