"""Filesystem helpers — chiefly, defeat Windows' 260-char MAX_PATH limit.

Deeply-nested NAS output trees (long experiment + gonad names mirrored several levels under a long
UNC base, e.g. ``\\\\NAS\\share\\...\\<exp>\\<exp>_nd2\\<gonad>\\<gonad>__image_summary.csv``) push the
longest filenames past 260 chars; pandas/tifffile then raise a *misleading* FileNotFoundError on a
directory that plainly exists. Wrapping the final path in the Windows extended-length form
(``\\\\?\\`` / ``\\\\?\\UNC\\``) lifts the limit. No-op off Windows.
"""
from __future__ import annotations

import os


def long_path(path) -> str:
    """Return a string path safe to pass to open()/pandas/tifffile past the 260-char MAX_PATH limit.

    On Windows, returns the extended-length form (``\\\\?\\C:\\...`` for local paths, ``\\\\?\\UNC\\server\\share\\...``
    for UNC/NAS paths); elsewhere — or if already prefixed — returns the path unchanged. The result is
    absolute and backslash-normalized, as the ``\\\\?\\`` prefix requires (it disables Win32 path parsing,
    so forward slashes or relative segments would otherwise break)."""
    s = os.fspath(path)
    if os.name != "nt" or s.startswith("\\\\?\\"):
        return s
    s = os.path.abspath(s)                       # \\?\ needs an absolute, separator-normalized path
    if s.startswith("\\\\"):                     # UNC: \\server\share\... -> \\?\UNC\server\share\...
        return "\\\\?\\UNC\\" + s[2:]
    return "\\\\?\\" + s
