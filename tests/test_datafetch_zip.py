"""germquant.datafetch, one member of a remote zip (design 6.1 and 6.9; package P3, acceptance 4): the end records
(classic and ZIP64), the central directory and one member are read with HTTP Range requests from a local server;
stored, deflated and bzip2 members come out whole with their CRC-32 checked; size, sha256 and CRC mismatches are
refused; a broken transfer resumes, stored or compressed; a member compressed with a method the reader cannot
decode (deflate64, LZMA) or encrypted is refused with the method named; a server that ignores Range is refused
rather than streamed. The archive is not trusted to describe itself: a pinned zip_size the server contradicts is
named, sizes that cannot be right are refused before any member data is requested, a member that inflates past its
recorded size stops one byte past it, and an oversized central directory is refused; bytes after the end record
are accepted, as Python's zipfile accepts them. Nothing leaves the machine."""
from __future__ import annotations

import ast
import bz2
import hashlib
import http.server
import io
import random
import re
import struct
import threading
import zipfile
import zlib
from contextlib import contextmanager
from pathlib import Path

import pytest

from germquant import datafetch as D

BIG = bytes((i * 7 + i // 251) % 256 for i in range(200_000))          # a member larger than the archive tail
TEXT = b"germquant zip member test\n" * 400
SMALL = b"a small member"
NOISE = random.Random(7).randbytes(150_000)                            # incompressible: deflate keeps its size


def build_zip(members, *, zip64: bool, encrypt: str | None = None) -> bytes:
    """A zip archive written by hand: ``members`` is [(name, data, method)]. With ``zip64`` every local and
    central entry carries its sizes (and the offset) in a ZIP64 extra field with 0xFFFFFFFF in the fixed fields,
    and the archive ends with a ZIP64 end record, its locator and an end record whose fields are saturated: the
    layout of an archive above 4 GB, at test size. A member named in ``encrypt`` gets the encryption flag."""
    out, central = bytearray(), bytearray()
    for name, data, method in members:
        nb = name.encode("utf-8")
        crc = zlib.crc32(data) & 0xFFFFFFFF
        if method == 8:
            c = zlib.compressobj(6, zlib.DEFLATED, -15)
            comp = c.compress(data) + c.flush()
        elif method == 12:
            comp = bz2.compress(data)
        else:
            comp = data                     # stored, or a method the reader refuses before it decodes anything
        flags = 0x800 | (0x1 if name == encrypt else 0)
        offset = len(out)
        if zip64:
            lextra = struct.pack("<HHQQ", 1, 16, len(data), len(comp))
            lsizes = (0xFFFFFFFF, 0xFFFFFFFF)
            cextra = struct.pack("<HHQQQ", 1, 24, len(data), len(comp), offset)
            csizes, coff = (0xFFFFFFFF, 0xFFFFFFFF), 0xFFFFFFFF
        else:
            lextra, lsizes, cextra, csizes, coff = b"", (len(comp), len(data)), b"", (len(comp), len(data)), offset
        out += struct.pack("<4sHHHHHIIIHH", b"PK\x03\x04", 45, flags, method, 0, 0x21, crc, lsizes[0], lsizes[1],
                           len(nb), len(lextra)) + nb + lextra + comp
        central += struct.pack("<4sHHHHHHIIIHHHHHII", b"PK\x01\x02", 45, 45, flags, method, 0, 0x21, crc, csizes[0],
                               csizes[1], len(nb), len(cextra), 0, 0, 0, 0, coff) + nb + cextra
    n, cd_off, cd_size = len(members), len(out), len(central)
    out += central
    if zip64:
        z64 = len(out)
        out += struct.pack("<4sQHHIIQQQQ", b"PK\x06\x06", 44, 45, 45, 0, 0, n, n, cd_size, cd_off)
        out += struct.pack("<4sIQI", b"PK\x06\x07", 0, z64, 1)
        out += struct.pack("<4sHHHHIIH", b"PK\x05\x06", 0, 0, 0xFFFF, 0xFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0)
    else:
        out += struct.pack("<4sHHHHIIH", b"PK\x05\x06", 0, 0, n, n, cd_size, cd_off, 0)
    return bytes(out)


class _RangeHandler(http.server.BaseHTTPRequestHandler):
    payload = b""
    ignore_range = False
    cut_at_start: int | None = None          # cut the first response for a range starting here ...
    cut_after = 0                            # ... after this many bytes
    log: list = []

    def log_message(self, *a):
        pass

    def do_GET(self):  # noqa: N802 - the http.server API
        cls, data = type(self), type(self).payload
        rng = self.headers.get("Range")
        cls.log.append(rng)
        if not rng or cls.ignore_range:
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        m = re.fullmatch(r"bytes=(\d*)-(\d*)", rng.strip())
        a, b = m.group(1), m.group(2)
        if a == "":
            start, end = len(data) - int(b), len(data) - 1
        else:
            start, end = int(a), (int(b) if b else len(data) - 1)
        if start >= len(data):                              # a range past the end: 416 with the file's size
            self.send_response(416)
            self.send_header("Content-Range", f"bytes */{len(data)}")
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        end = min(end, len(data) - 1)
        body = data[start:end + 1]
        self.send_response(206)
        self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if cls.cut_at_start is not None and start == cls.cut_at_start:
            cls.cut_at_start = None
            self.wfile.write(body[:cls.cut_after])
            self.wfile.flush()
            self.close_connection = True
            return
        self.wfile.write(body)


@contextmanager
def serve(payload: bytes, **opts):
    handler = type("H", (_RangeHandler,), {"payload": payload, "log": [], **opts})
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    try:
        yield f"http://127.0.0.1:{srv.server_address[1]}/archive.zip", handler
    finally:
        srv.shutdown()
        srv.server_close()


MEMBERS = [("folder/", b"", 0), ("folder/big.tif", BIG, 0), ("folder/text.txt", TEXT, 8),
           ("folder/small.bin", SMALL, 12)]


@pytest.mark.parametrize("zip64", [True, False])
def test_the_hand_built_archives_are_valid_zips(zip64):
    """The fixture itself: Python's zipfile reads both layouts (so the reader is tested on real zip files)."""
    blob = build_zip(MEMBERS, zip64=zip64)
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        assert zf.read("folder/big.tif") == BIG and zf.read("folder/text.txt") == TEXT
        assert zf.read("folder/small.bin") == SMALL
    assert (b"PK\x06\x06" in blob[-200:]) == zip64


@pytest.mark.parametrize("zip64", [True, False])
def test_index_and_members_come_out_whole(tmp_path, zip64):
    blob = build_zip(MEMBERS, zip64=zip64)
    with serve(blob) as (url, h):
        assert D.remote_size(url) == len(blob)
        idx = D.zip_index(url)
        assert sorted(idx) == sorted(m[0] for m in MEMBERS)
        big = idx["folder/big.tif"]
        assert (big.method_name, big.file_size, big.compressed_size) == ("stored", len(BIG), len(BIG))
        assert big.crc32 == zlib.crc32(BIG) & 0xFFFFFFFF and idx["folder/text.txt"].method_name == "deflate"
        for name, data, _ in MEMBERS[1:]:
            out = D.fetch_zip_member(url, name, tmp_path / name.split("/")[-1], size=len(data),
                                     sha256=hashlib.sha256(data).hexdigest(), crc32=f"{zlib.crc32(data) & 0xFFFFFFFF:08x}",
                                     zip_size=len(blob), index=idx)
            assert out.read_bytes() == data
            assert not D.part_path(out).exists() and not out.with_name(out.name + ".zpart").exists()
        assert all(r is not None for r in h.log)                         # every request asked for a range
        served = sum(int(e) - int(s) + 1 for s, e in (re.fullmatch(r"bytes=(\d+)-(\d+)", r).groups()
                                                     for r in h.log if re.fullmatch(r"bytes=\d+-\d+", r)))
        assert served < 2 * len(blob)                                   # never the whole archive per member


def test_a_large_archive_reads_its_central_directory_with_one_more_request(tmp_path):
    """An archive whose central directory lies before the tail that is read first (here: a member of about 1 MB
    behind the directory's start) takes one more Range request for the directory, as DeepFoci's 10 GB zip would
    if its directory were large."""
    filler = bytes(1_100_000)
    blob = build_zip([("a.bin", filler, 0)] + [(f"m{i:04d}.txt", b"x" * 10, 0) for i in range(3000)], zip64=True)
    with serve(blob) as (url, h):
        idx = D.zip_index(url, zip_size=len(blob))
        assert len(idx) == 3001
        out = D.fetch_zip_member(url, "m2999.txt", tmp_path / "m.txt", index=idx)
    assert out.read_bytes() == b"x" * 10


def test_size_sha256_and_crc_mismatches_are_refused(tmp_path):
    blob = build_zip(MEMBERS, zip64=True)
    with serve(blob) as (url, h):
        with pytest.raises(D.DownloadError, match="expected 12"):
            D.fetch_zip_member(url, "folder/text.txt", tmp_path / "a.txt", size=12)
        with pytest.raises(D.DownloadError, match="CRC-32"):
            D.fetch_zip_member(url, "folder/text.txt", tmp_path / "b.txt", crc32="00000000")
        assert not (tmp_path / "a.txt").exists() and not (tmp_path / "b.txt").exists()     # nothing downloaded
        with pytest.raises(D.DownloadError, match="sha256"):
            D.fetch_zip_member(url, "folder/big.tif", tmp_path / "c.tif", sha256="0" * 64)
        assert not (tmp_path / "c.tif").exists() and not D.part_path(tmp_path / "c.tif").exists()   # discarded
        with pytest.raises(D.ZipMemberError, match="no member 'nope.tif'.*folder/big.tif"):
            D.fetch_zip_member(url, "nope.tif", tmp_path / "d.tif")


@pytest.mark.parametrize("name, payload, method", [("folder/big.tif", BIG, 0), ("folder/noise.bin", NOISE, 8)],
                         ids=["stored", "deflated"])
def test_a_broken_transfer_resumes_with_a_range_request(tmp_path, name, payload, method):
    """Cut after 50,000 bytes of the member's data: the partial file (<dest>.part for a stored member, <dest>.zpart
    for the compressed bytes of a deflated one) stays, and the next call asks only for the rest."""
    blob = build_zip(MEMBERS + [("folder/noise.bin", NOISE, 8)], zip64=True)
    with serve(blob) as (url, h):
        idx = D.zip_index(url)
        m = idx[name]
        assert m.method == method and m.compressed_size > 50_000
        data_start = m.header_offset + 30 + len(name.encode()) + 20          # 20: the local ZIP64 extra field
        h.cut_at_start, h.cut_after = data_start, 50_000
        dest = tmp_path / name.split("/")[-1]
        with pytest.raises(D.DownloadError, match="run again to resume"):
            D.fetch_zip_member(url, name, dest, index=idx)
        part, zpart = D.part_path(dest), dest.with_name(dest.name + ".zpart")
        kept = part if method == 0 else zpart
        assert kept.stat().st_size == 50_000 and not dest.exists()
        assert not (zpart if method == 0 else part).exists()
        out = D.fetch_zip_member(url, name, dest, index=idx, sha256=hashlib.sha256(payload).hexdigest())
        assert h.log[-1] == f"bytes={data_start + 50_000}-{data_start + m.compressed_size - 1}"   # only the rest
    assert out.read_bytes() == payload and not part.exists() and not zpart.exists()


def patch_central(blob: bytes, name: str, *, csize: int | None = None, usize: int | None = None) -> bytes:
    """``blob`` (a classic archive from build_zip) with the central directory entry of ``name`` claiming other
    compressed and uncompressed sizes: an archive that lies about a member."""
    out = bytearray(blob)
    pos = out.find(b"PK\x01\x02")
    while pos >= 0:
        nlen = struct.unpack_from("<H", out, pos + 28)[0]
        if out[pos + 46:pos + 46 + nlen].decode() == name:
            if csize is not None:
                struct.pack_into("<I", out, pos + 20, csize)
            if usize is not None:
                struct.pack_into("<I", out, pos + 24, usize)
            return bytes(out)
        pos = out.find(b"PK\x01\x02", pos + 4)
    raise KeyError(name)


def test_member_sizes_that_cannot_be_right_are_refused_before_any_data(tmp_path):
    """A stored member whose two sizes differ, a deflated member far larger than its output and a member whose data
    would run past the end of the archive are refused from the central directory alone: no local header and no data
    is requested."""
    blob = build_zip([("a/stored.bin", SMALL, 0), ("a/text.txt", TEXT, 8)], zip64=False)
    cases = [
        (patch_central(blob, "a/stored.bin", csize=8 << 20), "a/stored.bin", "stored without compression"),
        (patch_central(blob, "a/text.txt", csize=5_000_000), "a/text.txt", "more than deflate ever produces"),
        (patch_central(blob, "a/stored.bin", csize=10_000_000, usize=10_000_000), "a/stored.bin",
         "would run past the end"),
    ]
    for lying, member, needle in cases:
        with serve(lying) as (url, h):
            idx = D.zip_index(url, zip_size=len(lying))
            n = len(h.log)
            with pytest.raises(D.ZipMemberError, match=needle):
                D.fetch_zip_member(url, member, tmp_path / "x.bin", index=idx, zip_size=len(lying))
            assert len(h.log) == n                                              # nothing more was requested
        assert not (tmp_path / "x.bin").exists() and not D.part_path(tmp_path / "x.bin").exists()


def test_a_member_that_inflates_past_its_recorded_size_stops_one_byte_past_it(tmp_path):
    """A deflate bomb: 8 MB of zeros (about 8 kB compressed) recorded as a 1,000-byte file. Inflating stops after
    1,001 bytes, and nothing is kept."""
    zeros = bytes(8 << 20)
    blob = patch_central(build_zip([("bomb.bin", zeros, 8)], zip64=False), "bomb.bin", usize=1000)
    with serve(blob) as (url, h):
        with pytest.raises(D.ZipMemberError, match=r"inflates past the 1000 bytes .*stopped after 1001 bytes"):
            D.fetch_zip_member(url, "bomb.bin", tmp_path / "bomb.bin")
    dest = tmp_path / "bomb.bin"
    assert not dest.exists() and not D.part_path(dest).exists() and not dest.with_name("bomb.bin.zpart").exists()
    src = tmp_path / "bomb.z"                                       # the inflater itself, on a file
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    src.write_bytes(c.compress(zeros) + c.flush())
    with pytest.raises(D.ZipMemberError, match="stopped after 1001 bytes"):
        D._decompress(src, tmp_path / "out.bin", 8, 1 << 16, 1000)
    src.write_bytes(bz2.compress(zeros))
    with pytest.raises(D.ZipMemberError, match="stopped after 1001 bytes"):
        D._decompress(src, tmp_path / "out.bin", 12, 1 << 16, 1000)
    assert not (tmp_path / "out.bin").exists()


@pytest.mark.parametrize("delta", [-1000, 1000, 10_000_000])
def test_a_pinned_zip_size_the_server_contradicts_is_named(delta):
    """A pinned zip_size that no longer matches the archive (a re-uploaded deposit) names both sizes, whether the pin
    is smaller, a little larger or far past the end (the server then answers 416)."""
    blob = build_zip(MEMBERS, zip64=True)
    with serve(blob) as (url, h):
        with pytest.raises(D.ZipMemberError, match=rf"reports {len(blob)} bytes, not the {len(blob) + delta} bytes"):
            D.zip_index(url, zip_size=len(blob) + delta)


@pytest.mark.parametrize("zip64", [True, False])
def test_bytes_after_the_end_record_are_accepted(tmp_path, zip64):
    blob = build_zip(MEMBERS, zip64=zip64) + b"\x00" * 16
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:                       # Python's zipfile reads it too
        assert zf.read("folder/text.txt") == TEXT
    with serve(blob) as (url, h):
        assert D.fetch_zip_member(url, "folder/text.txt", tmp_path / "t.txt").read_bytes() == TEXT


def test_an_oversized_central_directory_is_refused(monkeypatch):
    blob = build_zip(MEMBERS, zip64=True)
    monkeypatch.setattr(D, "MAX_CENTRAL_DIRECTORY", 100)
    with serve(blob) as (url, h):
        with pytest.raises(D.ZipMemberError, match="more than the 100 this reader holds in memory"):
            D.zip_index(url)


def test_bz2_is_imported_only_to_inflate_a_bzip2_member():
    """A Python built without _bz2 still imports datafetch, which the pack model path and pack status use."""
    tree = ast.parse(Path(D.__file__).read_text(encoding="utf-8"))
    top = {a.name for node in tree.body if isinstance(node, ast.Import) for a in node.names}
    top |= {node.module for node in tree.body if isinstance(node, ast.ImportFrom) and node.module}
    assert "bz2" not in top


@pytest.mark.parametrize("method, name", [(9, "deflate64"), (14, "LZMA"), (93, "zstd")])
def test_a_compression_it_cannot_read_is_named(tmp_path, method, name):
    blob = build_zip([("x/data.tif", BIG[:5000], method), ("x/ok.txt", TEXT, 8)], zip64=True)
    with serve(blob) as (url, h):
        with pytest.raises(D.UnsupportedZipCompression, match=rf"compressed with {name} \(zip method {method}\)"):
            D.fetch_zip_member(url, "x/data.tif", tmp_path / "data.tif")
        assert D.fetch_zip_member(url, "x/ok.txt", tmp_path / "ok.txt").read_bytes() == TEXT
    assert not (tmp_path / "data.tif").exists()


def test_an_encrypted_member_is_refused(tmp_path):
    blob = build_zip([("secret.bin", SMALL, 0)], zip64=False, encrypt="secret.bin")
    with serve(blob) as (url, h):
        with pytest.raises(D.UnsupportedZipCompression, match="encrypted"):
            D.fetch_zip_member(url, "secret.bin", tmp_path / "s.bin")


def test_a_server_that_ignores_range_is_refused(tmp_path):
    blob = build_zip(MEMBERS, zip64=True)
    with serve(blob, ignore_range=True) as (url, h):
        with pytest.raises(D.RangeNotSupported, match="Range"):
            D.zip_index(url)
        with pytest.raises(D.RangeNotSupported):
            D.fetch_zip_member(url, "folder/big.tif", tmp_path / "big.tif", zip_size=len(blob))


def test_not_a_zip_is_named():
    with serve(b"plain bytes, no end record" * 10) as (url, h):
        with pytest.raises(D.ZipMemberError, match="not a zip archive"):
            D.zip_index(url)


def test_python_zipfile_archives_with_zip64_entries_are_read(tmp_path):
    """An archive written by Python's zipfile with ZIP64 entries (force_zip64) and a deflated member."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        with zf.open("a/big.tif", "w", force_zip64=True) as fh:
            fh.write(BIG)
        zf.writestr(zipfile.ZipInfo("a/text.txt"), TEXT, compress_type=zipfile.ZIP_DEFLATED)
    blob = buf.getvalue()
    with serve(blob) as (url, h):
        assert D.fetch_zip_member(url, "a/big.tif", tmp_path / "big.tif", size=len(BIG)).read_bytes() == BIG
        assert D.fetch_zip_member(url, "a/text.txt", tmp_path / "t.txt").read_bytes() == TEXT
