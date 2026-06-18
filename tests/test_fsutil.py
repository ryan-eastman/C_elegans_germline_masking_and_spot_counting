"""long_path: defeat the Windows 260-char MAX_PATH limit that deeply-nested NAS output trees hit
(the bug that failed ..._H001__image_summary.csv mid-batch)."""
import os

import pytest

from germquant.fsutil import long_path


@pytest.mark.skipif(os.name != "nt", reason="extended-length prefixes are Windows-only")
def test_long_path_prefixes():
    assert long_path(r"C:\a\b").startswith("\\\\?\\C:\\")
    assert long_path(r"\\NAS\share\a").startswith("\\\\?\\UNC\\NAS\\share")
    assert long_path("\\\\?\\C:\\already") == "\\\\?\\C:\\already"   # idempotent


def test_long_path_noop_off_windows():
    if os.name == "nt":
        pytest.skip("posix-only check")
    assert long_path("/tmp/a/b") == "/tmp/a/b"


def test_write_path_over_260_chars(tmp_path):
    r"""Build a path well past 260 chars and confirm mkdir + write + read all succeed through
    long_path (on Windows this exercises the \\?\ mechanism; elsewhere it's just a deep path)."""
    deep = tmp_path
    for _ in range(6):
        deep = deep / ("g" * 40)
    target = deep / ("f" * 40 + "__image_summary.csv")
    assert len(str(target)) > 260

    os.makedirs(long_path(deep), exist_ok=True)
    with open(long_path(target), "w", encoding="utf-8") as fh:
        fh.write("ok")
    with open(long_path(target), "r", encoding="utf-8") as fh:
        assert fh.read() == "ok"
