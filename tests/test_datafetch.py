"""germquant.datafetch: resumable downloads with size and sha256 checks, against a local HTTP server that
honours Range requests (and one that cuts the first transfer short). Nothing leaves the machine."""
from __future__ import annotations

import hashlib
import http.server
import threading
from contextlib import contextmanager

import pytest

from germquant import datafetch as D

PAYLOAD = bytes(range(256)) * 400           # 102,400 bytes


class _Handler(http.server.BaseHTTPRequestHandler):
    payload = PAYLOAD
    cut_first_at: int | None = None          # send this many bytes of the first full GET, then drop the connection
    ignore_range = False
    range_from_zero = False                  # answer a Range request with 206 but the whole file (bytes 0-...)
    log: list = []

    def log_message(self, *a):                # quiet
        pass

    def do_GET(self):  # noqa: N802 - the http.server API
        data = type(self).payload
        rng = self.headers.get("Range")
        type(self).log.append((self.path, rng))
        start = 0
        if rng and not type(self).ignore_range:
            start = 0 if type(self).range_from_zero else int(rng.split("=", 1)[1].split("-", 1)[0])
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
        else:
            self.send_response(200)
        body = data[start:]
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        cut = type(self).cut_first_at
        if cut is not None and not rng:
            type(self).cut_first_at = None
            self.wfile.write(body[:cut])
            self.wfile.flush()
            self.close_connection = True
            return
        self.wfile.write(body)


@contextmanager
def serve(payload=PAYLOAD, cut_first_at=None, ignore_range=False, range_from_zero=False):
    handler = type("H", (_Handler,), {"payload": payload, "cut_first_at": cut_first_at,
                                      "ignore_range": ignore_range, "range_from_zero": range_from_zero, "log": []})
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    try:
        yield f"http://127.0.0.1:{srv.server_address[1]}/blob", handler
    finally:
        srv.shutdown()
        srv.server_close()


SHA = hashlib.sha256(PAYLOAD).hexdigest()


def test_download_checks_size_and_sha(tmp_path):
    with serve() as (url, h):
        out = D.download(url, tmp_path / "a" / "blob", size=len(PAYLOAD), sha256=SHA)
    assert out.read_bytes() == PAYLOAD and not D.part_path(out).exists()
    with serve() as (url, h):
        with pytest.raises(D.DownloadError, match="sha256"):
            D.download(url, tmp_path / "b", size=len(PAYLOAD), sha256="0" * 64)
    assert not (tmp_path / "b").exists() and not D.part_path(tmp_path / "b").exists()     # discarded
    with serve() as (url, h):
        with pytest.raises(D.DownloadError, match="ended early"):
            D.download(url, tmp_path / "c", size=len(PAYLOAD) + 10, sha256=SHA)
    assert not (tmp_path / "c").exists()


def test_download_resumes_with_a_range_request(tmp_path):
    dest = tmp_path / "blob"
    with serve(cut_first_at=30_000) as (url, h):
        with pytest.raises(Exception):                                    # the connection drops mid-way
            D.download(url, dest, size=len(PAYLOAD), sha256=SHA)
        part = D.part_path(dest)
        assert part.is_file() and 0 < part.stat().st_size < len(PAYLOAD) and not dest.exists()
        have = part.stat().st_size
        D.download(url, dest, size=len(PAYLOAD), sha256=SHA)
        assert h.log[-1] == ("/blob", f"bytes={have}-")                    # it asked for the rest only
    assert dest.read_bytes() == PAYLOAD and not part.exists()


def test_a_server_that_ignores_range_is_read_from_the_start(tmp_path):
    dest = tmp_path / "blob"
    D.part_path(dest).write_bytes(PAYLOAD[:1000])
    with serve(ignore_range=True) as (url, h):
        D.download(url, dest, size=len(PAYLOAD), sha256=SHA)
    assert dest.read_bytes() == PAYLOAD


def test_a_206_that_starts_elsewhere_is_read_from_the_start(tmp_path):
    """A server that answers the Range request with 206 but Content-Range bytes 0-...: appending would double
    the first bytes, so the download starts again from byte 0 and the file comes out whole."""
    dest = tmp_path / "blob"
    D.part_path(dest).write_bytes(PAYLOAD[:1000])
    with serve(range_from_zero=True) as (url, h):
        D.download(url, dest, size=len(PAYLOAD), sha256=SHA)
        assert h.log[-1] == ("/blob", "bytes=1000-")
    assert dest.read_bytes() == PAYLOAD and not D.part_path(dest).exists()
