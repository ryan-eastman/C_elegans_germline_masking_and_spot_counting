"""Resumable downloads with size and sha256 checks (docs/DESIGN_system_packs.md, sections 3.5 and 6.9).

``download(url, dest, size=..., sha256=...)`` streams to ``<dest>.part`` and asks the server to continue a
partial ``.part`` file with an HTTP Range request (``Range: bytes=<n>-``); a server that ignores the range
(status 200 instead of 206) is read again from the start. Only a file whose size and sha256 match is moved
to ``dest`` (os.replace, so ``dest`` never holds a partial file); a mismatch deletes the ``.part`` file and
raises DownloadError, and a broken connection keeps it for the next call. Standard library only (urllib),
so it runs in every environment, CI included, and imports nothing heavy.

One member of a remote zip (package P3). Public datasets such as DeepFoci ship as one zip of 10 to 50 GB;
``fetch_zip_member(url, member, dest, ...)`` reads only what one member needs, with HTTP Range requests: the
end of the archive (the end-of-central-directory record, and the ZIP64 locator and record of an archive
above 4 GB or 65535 entries), the central directory, the member's local header and its data. Members stored
without compression (method 0), deflated (8) or bzip2-compressed (12) are read; any other method (deflate64,
LZMA, zstd and the rest) or an encrypted member raises UnsupportedZipCompression, which names the method.
The member's CRC-32 from the central directory is always checked, and ``size`` and ``sha256`` when given;
a partial transfer is kept and resumed like a plain download. A server that does not honour Range requests
is refused (RangeNotSupported) instead of streaming the whole archive.

The pack model fetcher (``germquant pack fetch-models``) and the example fetcher (``germquant pack
fetch-example``, germquant.packs.examples) use it.
"""
from __future__ import annotations

import bz2
import hashlib
import logging
import os
import struct
import time
import urllib.request
import zlib
from collections.abc import Callable
from dataclasses import dataclass
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


# ---- HTTP Range reads, and one member of a remote zip ---------------------------------------------------------

class RangeNotSupported(RuntimeError):
    """A server that answers a Range request with the whole file or with another range: reading one member of a
    large archive then means streaming the whole archive, which this reader refuses to do."""


class ZipMemberError(RuntimeError):
    """A remote zip whose end records or central directory cannot be read, or that has no such member."""


class UnsupportedZipCompression(ZipMemberError):
    """A zip member compressed with a method this reader cannot decode, or an encrypted member. The message
    names the method (deflate64, LZMA, zstd, ...)."""


# the compression methods of the zip specification (APPNOTE 4.4.5), by number, for messages
ZIP_METHODS = {0: "stored", 1: "shrink", 2: "reduce (factor 1)", 3: "reduce (factor 2)", 4: "reduce (factor 3)",
               5: "reduce (factor 4)", 6: "implode", 8: "deflate", 9: "deflate64", 10: "PKWARE DCL implode",
               12: "bzip2", 14: "LZMA", 16: "IBM z/OS CMPSC", 18: "IBM TERSE", 19: "IBM LZ77",
               20: "zstd (old id)", 93: "zstd", 94: "MP3", 95: "xz", 96: "JPEG", 97: "WavPack", 98: "PPMd",
               99: "AES encryption"}
READABLE_METHODS = (0, 8, 12)              # stored, deflate (zlib) and bzip2 (bz2): the standard library's
_EOCD_SIG, _LOC64_SIG, _EOCD64_SIG = b"PK\x05\x06", b"PK\x06\x07", b"PK\x06\x06"
_CDIR_SIG, _LOCAL_SIG = b"PK\x01\x02", b"PK\x03\x04"
_EOCD_LEN, _LOC64_LEN, _EOCD64_LEN, _CDIR_LEN, _LOCAL_LEN = 22, 20, 56, 46, 30
_TAIL = _EOCD_LEN + 0xFFFF + _LOC64_LEN + _EOCD64_LEN   # the end record with the longest comment, plus ZIP64
_U16, _U32 = 0xFFFF, 0xFFFFFFFF


@dataclass(frozen=True)
class ZipMember:
    """One entry of a zip's central directory: sizes and the local header offset, ZIP64 values applied."""

    name: str
    method: int
    flags: int
    compressed_size: int
    file_size: int
    crc32: int
    header_offset: int

    @property
    def method_name(self) -> str:
        return ZIP_METHODS.get(self.method, f"unknown method {self.method}")

    @property
    def encrypted(self) -> bool:
        return bool(self.flags & 0x1)

    @property
    def is_dir(self) -> bool:
        return self.name.endswith("/")


def _open(url: str, headers: dict | None, timeout: float):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
    return urllib.request.urlopen(req, timeout=timeout)


def _content_range(header: str | None) -> tuple[int, int, int | None] | None:
    """(first, last, total) of a ``Content-Range: bytes a-b/n`` header (total None for ``*``), or None."""
    if not header:
        return None
    try:
        unit, rest = header.strip().split(" ", 1)
        if unit.lower() != "bytes":
            return None
        span, total = rest.split("/", 1)
        first, last = span.split("-", 1)
        return int(first), int(last), (None if total.strip() == "*" else int(total))
    except (ValueError, IndexError):
        return None


def _no_range(url: str, status, what: str) -> RangeNotSupported:
    return RangeNotSupported(f"{url}: the server did not answer the HTTP Range request for {what} with that range "
                             f"(status {status}); one member of a remote archive can be read only from a server "
                             f"that honours Range requests (download the whole archive instead)")


def remote_size(url: str, *, timeout: float = 120.0, headers: dict | None = None) -> int:
    """The size of a remote file, from the Content-Range of a one-byte Range request (which also shows that the
    server honours ranges). RangeNotSupported otherwise."""
    with _open(url, {**(headers or {}), "Range": "bytes=0-0"}, timeout) as r:
        status = getattr(r, "status", None) or r.getcode()
        cr = _content_range(r.headers.get("Content-Range"))
        if status != 206 or cr is None or cr[0] != 0 or cr[2] is None:
            raise _no_range(url, status, "its first byte")
        r.read(1)
        return int(cr[2])


def read_range(url: str, start: int, length: int, *, timeout: float = 120.0, headers: dict | None = None) -> bytes:
    """Bytes ``[start, start + length)`` of ``url`` through one Range request. RangeNotSupported when the
    server answers with another range or the whole file; ZipMemberError when fewer bytes arrive."""
    if length <= 0:
        return b""
    with _open(url, {**(headers or {}), "Range": f"bytes={start}-{start + length - 1}"}, timeout) as r:
        status = getattr(r, "status", None) or r.getcode()
        cr = _content_range(r.headers.get("Content-Range"))
        if status != 206 or cr is None or cr[0] != start:
            raise _no_range(url, status, f"bytes {start} to {start + length - 1}")
        data = r.read(length)
    if len(data) != length:
        raise ZipMemberError(f"{url}: {len(data)} of the {length} bytes at offset {start} arrived")
    return data


def _find_eocd(tail: bytes) -> int:
    """The offset in ``tail`` (the end of the archive) of the end-of-central-directory record: the last
    signature whose comment length reaches exactly the end of the file."""
    end = len(tail) - _EOCD_LEN + 4
    while end >= 4:
        i = tail.rfind(_EOCD_SIG, 0, end)
        if i < 0:
            break
        if i + _EOCD_LEN <= len(tail):
            (comment_len,) = struct.unpack_from("<H", tail, i + 20)
            if i + _EOCD_LEN + comment_len == len(tail):
                return i
        end = i + 3
    raise ZipMemberError("no end-of-central-directory record: not a zip archive, or a truncated one")


def _zip64_extra(extra: bytes, usize: int, csize: int, offset: int, disk: int) -> tuple[int, int, int, int]:
    """Apply the ZIP64 extended information field (header 0x0001) of a central directory entry: each value
    saturated in the fixed fields (0xFFFFFFFF, or 0xFFFF for the disk) is read from it, in the order
    uncompressed size, compressed size, local header offset, disk."""
    pos = 0
    while pos + 4 <= len(extra):
        head, n = struct.unpack_from("<HH", extra, pos)
        body = extra[pos + 4:pos + 4 + n]
        if head == 0x0001:
            q = 0
            values = []
            for saturated, width in ((usize == _U32, 8), (csize == _U32, 8), (offset == _U32, 8), (disk == _U16, 4)):
                if not saturated:
                    values.append(None)
                    continue
                if q + width > len(body):
                    raise ZipMemberError("a ZIP64 extra field is shorter than the values it must hold")
                values.append(struct.unpack_from("<Q" if width == 8 else "<I", body, q)[0])
                q += width
            usize = values[0] if values[0] is not None else usize
            csize = values[1] if values[1] is not None else csize
            offset = values[2] if values[2] is not None else offset
            disk = values[3] if values[3] is not None else disk
            break
        pos += 4 + n
    return usize, csize, offset, disk


def parse_central_directory(cd: bytes) -> list[ZipMember]:
    """The entries of a zip central directory (the bytes between its offset and offset + size)."""
    out: list[ZipMember] = []
    pos = 0
    while pos + _CDIR_LEN <= len(cd):
        if cd[pos:pos + 4] != _CDIR_SIG:
            raise ZipMemberError(f"central directory entry {len(out)}: bad signature")
        (_sig, _made, _need, flags, method, _time, _date, crc, csize, usize, nlen, xlen, clen, disk, _iattr,
         _eattr, offset) = struct.unpack_from("<4sHHHHHHIIIHHHHHII", cd, pos)
        raw_name = cd[pos + _CDIR_LEN:pos + _CDIR_LEN + nlen]
        name = raw_name.decode("utf-8" if flags & 0x800 else "cp437", "replace")
        extra = cd[pos + _CDIR_LEN + nlen:pos + _CDIR_LEN + nlen + xlen]
        usize, csize, offset, disk = _zip64_extra(extra, usize, csize, offset, disk)
        if disk != 0:
            raise ZipMemberError(f"member {name!r} starts on disk {disk}: spanned archives are not supported")
        out.append(ZipMember(name=name, method=method, flags=flags, compressed_size=int(csize), file_size=int(usize),
                             crc32=int(crc), header_offset=int(offset)))
        pos += _CDIR_LEN + nlen + xlen + clen
    return out


def _slice_or_read(url: str, tail: bytes, tail_start: int, offset: int, length: int, **kw) -> bytes:
    """Bytes ``[offset, offset + length)`` of the archive: from the tail already read when they lie inside it,
    else through one more Range request."""
    if offset >= tail_start and offset + length <= tail_start + len(tail):
        return tail[offset - tail_start:offset - tail_start + length]
    return read_range(url, offset, length, **kw)


def zip_index(url: str, *, zip_size: int | None = None, timeout: float = 120.0,
              headers: dict | None = None) -> dict[str, ZipMember]:
    """``{name: ZipMember}`` of the remote zip at ``url``, read with Range requests: the end of the archive,
    the ZIP64 record when the locator points at one, and the central directory. ``zip_size`` (the archive's
    size, when known) saves the request that asks for it."""
    kw = {"timeout": timeout, "headers": headers}
    total = int(zip_size) if zip_size else remote_size(url, **kw)
    n = min(total, _TAIL)
    tail_start = total - n
    tail = read_range(url, tail_start, n, **kw)
    i = _find_eocd(tail)
    disk, cd_disk, _n_disk, n_entries, cd_size, cd_offset = struct.unpack_from("<HHHHII", tail, i + 4)
    j = i - _LOC64_LEN
    if j >= 0 and tail[j:j + 4] == _LOC64_SIG:
        z64_disk, z64_offset, n_disks = struct.unpack_from("<IQI", tail, j + 4)
        if z64_disk != 0 or n_disks > 1:
            raise ZipMemberError(f"{url}: a spanned (multi-disk) archive is not supported")
        rec = _slice_or_read(url, tail, tail_start, z64_offset, _EOCD64_LEN, **kw)
        if rec[:4] != _EOCD64_SIG:
            raise ZipMemberError(f"{url}: the ZIP64 locator points at offset {z64_offset}, where no ZIP64 end record is")
        (_sig, _rec_size, _made, _need, disk, cd_disk, _n_disk, n_entries, cd_size,
         cd_offset) = struct.unpack("<4sQHHIIQQQQ", rec)
    elif n_entries == _U16 or cd_size == _U32 or cd_offset == _U32:
        raise ZipMemberError(f"{url}: the end record points to ZIP64 values, but no ZIP64 locator precedes it")
    if disk != 0 or cd_disk != 0:
        raise ZipMemberError(f"{url}: a spanned (multi-disk) archive is not supported")
    if cd_offset + cd_size > total:
        raise ZipMemberError(f"{url}: the central directory ({cd_size} bytes at {cd_offset}) lies past the end of "
                             f"the {total}-byte archive")
    members = parse_central_directory(_slice_or_read(url, tail, tail_start, cd_offset, cd_size, **kw))
    if len(members) != n_entries:
        raise ZipMemberError(f"{url}: the central directory holds {len(members)} entries, the end record says "
                             f"{n_entries}")
    return {m.name: m for m in members}


def _stream_range(url: str, part: Path, start: int, length: int, *, resume: bool, timeout: float, chunk: int,
                  progress: Callable[[int, int | None], None] | None, headers: dict | None) -> None:
    """Write bytes ``[start, start + length)`` of ``url`` to ``part``, continuing a partial ``part`` when
    ``resume``. A broken connection leaves ``part`` for the next call; RangeNotSupported when the server does
    not answer with the range asked for."""
    have = part.stat().st_size if (resume and part.is_file()) else 0
    if not resume and part.exists():
        part.unlink()
    if have > length:
        log.warning("%s: partial file larger than the expected %d bytes; starting again", part, length)
        part.unlink()
        have = 0
    if have == length:
        return
    first, last = start + have, start + length - 1
    t0, at_start = time.time(), have
    with _open(url, {**(headers or {}), "Range": f"bytes={first}-{last}"}, timeout) as r:
        status = getattr(r, "status", None) or r.getcode()
        cr = _content_range(r.headers.get("Content-Range"))
        if status != 206 or cr is None or cr[0] != first:
            raise _no_range(url, status, f"bytes {first} to {last}")
        with open(part, "ab" if have else "wb") as out:
            while have < length:
                block = r.read(min(chunk, length - have))
                if not block:
                    break
                out.write(block)
                have += len(block)
                if progress is not None:
                    progress(have, length)
    log.info("%s: %d bytes in %.1f s", url, have - at_start, time.time() - t0)
    if have < length:
        raise DownloadError(f"{url}: the transfer ended early ({have} of {length} bytes); run again to resume")


def _crc32_and_sha256(path: Path, chunk: int = 1 << 24) -> tuple[int, str]:
    crc, h = 0, hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            crc = zlib.crc32(block, crc)
            h.update(block)
    return crc & _U32, h.hexdigest()


def _decompress(src: Path, dest: Path, method: int, chunk: int) -> None:
    """Inflate a deflate (8) or bzip2 (12) stream from ``src`` into ``dest``; ZipMemberError when it is corrupt
    or ends early."""
    d = zlib.decompressobj(-15) if method == 8 else bz2.BZ2Decompressor()
    try:
        with open(src, "rb") as fin, open(dest, "wb") as out:
            for block in iter(lambda: fin.read(chunk), b""):
                out.write(d.decompress(block))
            if method == 8:
                out.write(d.flush())
    except (zlib.error, OSError, EOFError) as e:
        dest.unlink(missing_ok=True)
        raise ZipMemberError(f"{src.name}: the member's {ZIP_METHODS[method]} data is corrupt ({e})") from e
    if not d.eof:
        dest.unlink(missing_ok=True)
        raise ZipMemberError(f"{src.name}: the member's {ZIP_METHODS[method]} data ends early")


def fetch_zip_member(url: str, member: str, dest: str | Path, *, size: int | None = None, sha256: str | None = None,
                     crc32: str | int | None = None, zip_size: int | None = None,
                     index: dict[str, ZipMember] | None = None, resume: bool = True, timeout: float = 120.0,
                     chunk: int = CHUNK, progress: Callable[[int, int | None], None] | None = None,
                     headers: dict | None = None) -> Path:
    """Fetch the member ``member`` of the remote zip at ``url`` to ``dest`` with Range requests (module
    docstring). ``size`` (the member's uncompressed size), ``sha256`` and ``crc32`` (8 hex digits or an int)
    are checked when given; the member's CRC-32 from the central directory is always checked. ``zip_size``
    (the archive's size) and ``index`` (a zip_index already read) save requests. A transfer that breaks keeps
    ``<dest>.part`` (``<dest>.zpart`` for a compressed member) for the next call; a mismatch discards it and
    raises DownloadError. Returns ``dest``."""
    dest = Path(dest)
    kw = {"timeout": timeout, "headers": headers}
    idx = index if index is not None else zip_index(url, zip_size=zip_size, **kw)
    m = idx.get(member)
    if m is None or m.is_dir:
        files = sorted(n for n in idx if not n.endswith("/"))
        raise ZipMemberError(f"{url}: no member {member!r} (the archive holds {len(files)} files, for example "
                             f"{', '.join(files[:4]) or 'none'})")
    if m.encrypted:
        raise UnsupportedZipCompression(f"{url}: member {member!r} is encrypted, which this reader cannot read")
    if m.method not in READABLE_METHODS:
        raise UnsupportedZipCompression(
            f"{url}: member {member!r} is compressed with {m.method_name} (zip method {m.method}), which this reader "
            f"cannot decode; it reads stored, deflate and bzip2 members. Download the whole archive and extract it "
            f"with a tool that reads {m.method_name}, or take the example's documented fallback")
    if size is not None and int(size) != m.file_size:
        raise DownloadError(f"{url}: member {member!r} is {m.file_size} bytes in the archive, expected {size}; "
                            f"nothing was downloaded")
    if crc32 is not None:
        want_crc = int(crc32, 16) if isinstance(crc32, str) else int(crc32)
        if want_crc != m.crc32:
            raise DownloadError(f"{url}: member {member!r} has CRC-32 {m.crc32:08x} in the archive, expected "
                                f"{want_crc:08x}; nothing was downloaded")
    local = read_range(url, m.header_offset, _LOCAL_LEN, **kw)
    if local[:4] != _LOCAL_SIG:
        raise ZipMemberError(f"{url}: no local header for member {member!r} at offset {m.header_offset}")
    nlen, xlen = struct.unpack_from("<HH", local, 26)
    data_start = m.header_offset + _LOCAL_LEN + nlen + xlen
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = part_path(dest)
    stream = dict(resume=resume, timeout=timeout, chunk=chunk, progress=progress, headers=headers)
    if m.method == 0:
        _stream_range(url, part, data_start, m.compressed_size, **stream)
    else:
        zpart = dest.with_name(dest.name + ".zpart")
        _stream_range(url, zpart, data_start, m.compressed_size, **stream)
        try:
            _decompress(zpart, part, m.method, chunk)
        finally:
            zpart.unlink(missing_ok=True)
    got_size = part.stat().st_size
    if got_size != m.file_size:
        part.unlink(missing_ok=True)
        raise DownloadError(f"{url}: member {member!r} came out as {got_size} bytes, expected {m.file_size}; the "
                            f"file was discarded")
    got_crc, got_sha = _crc32_and_sha256(part)
    if got_crc != m.crc32:
        part.unlink(missing_ok=True)
        raise DownloadError(f"{url}: member {member!r} has CRC-32 {got_crc:08x}, the archive says {m.crc32:08x}; the "
                            f"file was discarded")
    if sha256 and got_sha != str(sha256).lower():
        part.unlink(missing_ok=True)
        raise DownloadError(f"{url}: member {member!r} has sha256 {got_sha[:12]}..., expected {str(sha256)[:12]}...; "
                            f"the file was discarded")
    os.replace(part, dest)
    return dest
