"""Spot quantification (RAD-51 etc.) via SpotMAX — replaces the v1 blob_log foci detector."""
from .detect import detect_spots

__all__ = ["detect_spots"]
