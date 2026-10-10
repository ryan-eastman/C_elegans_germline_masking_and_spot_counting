"""Spot quantification (RAD-51 etc.) through pluggable detection backends; SpotMAX is the validated
default (spots/backends.py, docs/BACKENDS.md)."""
from .backends import discover_backends, get_backend, list_backends, register_backend
from .detect import detect_spots
from .export import spots_to_image

__all__ = ["detect_spots", "discover_backends", "get_backend", "list_backends", "register_backend",
           "spots_to_image"]
