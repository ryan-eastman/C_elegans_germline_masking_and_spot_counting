"""Resumable downloads with size and sha256 checks (docs/DESIGN_system_packs.md, sections 3.5 and 6.9).

``download(url, dest, size=..., sha256=...)`` streams to ``<dest>.part`` and asks the server to continue a
partial ``.part`` file with an HTTP Range request (``Range: bytes=<n>-``); a server that ignores the range
(status 200 instead of 206) is read again from the start. Only a file whose size and sha256 match is moved
to ``dest`` (os.replace, so ``dest`` never holds a partial file); a mismatch deletes the ``.part`` file and
raises DownloadError, and a broken connection keeps it for the next call. Standard library only (urllib),
so it runs in every environment, CI included, and imports nothing heavy.

The pack model fetcher (``germquant pack fetch-models``) and, from package P3, the example fetcher use it.
"""
from __future__ import annotations

import hashlib
import logging
import os
import time
import urllib.request
from collections.abc import Callable
from pathlib import Path

log = logging.getLogger(__name__)

CHUNK = 1 << 22
USER_AGENT = "germquant datafetch"


class DownloadError(RuntimeError):
    """A download that ended with the wrong size or sha256 (the partial file is discarded)."""


def sha256_file(path: str | Path, chunk: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def part_path(dest: str | Path) -> Path:
    dest = Path(dest)
    return dest.with_name(dest.name + ".part")


def _content_range_start(header: str | None) -> int | None:
    """The first byte of a ``Content-Range: bytes a-b/n`` header, or None."""
    if not header:
        return None
    try:
        unit, rest = header.strip().split(" ", 1)
        if unit.lower() != "bytes":
            return None
        return int(rest.split("-", 1)[0])
    except (ValueError, IndexError):
        return None


def download(url: str, dest: str | Path, *, size: int | None = None, sha256: str | None = None,
             resume: bool = True, timeout: float = 120.0, chunk: int = CHUNK,
             progress: Callable[[int, int | None], None] | None = None, headers: dict | None = None) -> Path:
    """Fetch ``url`` to ``dest``, resuming ``<dest>.part`` when ``resume`` and the server honours a Range
    request, and verifying ``size`` and ``sha256`` (each when given) before ``dest`` appears. Returns ``dest``.

    Raises DownloadError when the finished file has the wrong size or digest (the partial file is deleted),
    and the urllib/OSError of a broken connection (the partial file is kept, so the next call resumes)."""
    dest = Path(dest)
    part = part_path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    have = part.stat().st_size if (resume and part.is_file()) else 0
    if not resume and part.exists():
        part.unlink()
    if size is not None and have > size:
        log.warning("%s: partial file larger than the expected %d bytes; starting again", part, size)
        part.unlink()
        have = 0
    if size is None or have < size:
        hdrs = {"User-Agent": USER_AGENT, **(headers or {})}
        if have:
            hdrs["Range"] = f"bytes={have}-"
        req = urllib.request.Request(url, headers=hdrs)
        t0, done_at_start = time.time(), have
        with urllib.request.urlopen(req, timeout=timeout) as r:
            status = getattr(r, "status", None) or r.getcode()
            if have and (status != 206 or _content_range_start(r.headers.get("Content-Range")) != have):
                log.info("%s: the server did not continue at byte %d (status %s); reading from the start",
                         url, have, status)
                have = 0
            mode = "ab" if have else "wb"
            with open(part, mode) as out:
                while True:
                    block = r.read(chunk)
                    if not block:
                        break
                    out.write(block)
                    have += len(block)
                    if size is not None and have > size:
                        break
                    if progress is not None:
                        progress(have, size)
        log.info("%s: %d bytes in %.1f s", url, have - done_at_start, time.time() - t0)
    got_size = part.stat().st_size if part.is_file() else 0
    if size is not None and got_size != size:
        if got_size > size:
            part.unlink(missing_ok=True)
            raise DownloadError(f"{url}: got {got_size} bytes, expected {size}; the file was discarded")
        raise DownloadError(f"{url}: the download ended early ({got_size} of {size} bytes); run again to resume")
    if sha256:
        got = sha256_file(part)
        if got != str(sha256).lower():
            part.unlink(missing_ok=True)
            raise DownloadError(f"{url}: sha256 {got[:12]}... does not match the expected {str(sha256)[:12]}...; "
                                f"the file was discarded")
    os.replace(part, dest)
    return dest
