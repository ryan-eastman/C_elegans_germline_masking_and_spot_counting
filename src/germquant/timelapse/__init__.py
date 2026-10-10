"""Time-lapse quantification (docs/DESIGN_timelapse.md, docs/TIMELAPSE.md).

The per-timepoint loop lives in ``series`` (process_series: every enabled stage on each timepoint of one
time-lapse file, the series record, the series_summary table and the reference volumes for registration);
the drift registration in ``register`` (register_series, ``germquant register``: the register_shifts table and
the registration report); ``layout`` reads a finished series folder back (timepoints, label images, nuclei
tables, the shift table, registered coordinates); ``track`` is the tracking exchange (the Cell-ACDC export, the
overlap / laptrack / none linking backends, the TrackMate / Imaris / Cell-ACDC importers, the track_points and
tracks tables, the strips) and ``track_validate`` the biologist validation (the stratified pick of hand-followed
nuclei and ``germquant track-validate``).
"""
from .layout import SeriesLayout, SeriesLayoutError, find_series, load_shifts, object_table
from .register import ReferenceVolumesMissing, register_series
from .series import SERIES_RECORD_SUFFIX, parse_timepoints, process_series, read_series_record
from .track import export_cellacdc, import_tracks, link_series, track_strips
from .track_validate import pick_validation_nuclei, validate_tracks

__all__ = [
    "SERIES_RECORD_SUFFIX", "ReferenceVolumesMissing", "SeriesLayout", "SeriesLayoutError", "export_cellacdc",
    "find_series", "import_tracks", "link_series", "load_shifts", "object_table", "parse_timepoints",
    "pick_validation_nuclei", "process_series", "read_series_record", "register_series", "track_strips",
    "validate_tracks",
]
