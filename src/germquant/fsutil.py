"""Filesystem helpers — chiefly, defeat Windows' 260-char MAX_PATH limit.

Deeply-nested NAS output trees (long experiment + gonad names mirrored several levels under a long
UNC base, e.g. ``\\\\NAS\\share\\...\\<exp>\\<exp>_nd2\\<gonad>\\<gonad>__image_summary.csv``) push the
longest filenames past 260 chars; pandas/tifffile then raise a *misleading* FileNotFoundError on a
directory that plainly exists. Wrapping the final path in the Windows extended-length form
(``\\\\?\\`` / ``\\\\?\\UNC\\``) lifts the limit. No-op off Windows.
"""
from __future__ import annotations

import os
from pathlib import Path

# The codec the config and channel-map files are decoded with, both for parsing and for the
# config_hash that goes into every manifest and table. The files are UTF-8 with a few non-ASCII
# characters in comments; before this constant existed they were read with the machine's locale
# codec, which on the validated Windows workstation is cp1252 (those characters decode to other
# characters there, harmless for YAML but part of the hashed text). Every golden, every analysis
# table and every run manifest carries hashes made that way, so the codec is pinned here rather
# than left to the locale: Linux, macOS and Python's UTF-8 mode now produce the workstation's hash
# instead of their own. Switching this to "utf-8" is the right long-term choice, but it changes every
# config_hash (and with it the resume bookkeeping of existing batch outputs), so do it only as a
# deliberate provenance step: regenerate tests/golden/*/config_hash.txt and note it in the roadmap.
CONFIG_TEXT_ENCODING = "cp1252"


def read_config_text(path) -> str:
    """The text of a config or channel-map file, decoded with ``CONFIG_TEXT_ENCODING`` on every
    platform (see the note above). Universal newlines, like ``Path.read_text``."""
    return Path(path).read_text(encoding=CONFIG_TEXT_ENCODING)


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
